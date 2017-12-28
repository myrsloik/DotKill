#include "VapourSynth.h"
#include "VSHelper.h"
#include <algorithm>
#include <memory>

typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;
    int napply;
    bool ignoreMatch;
} DotKillData;

template<typename T>
static T clamp(T v, T x, T y) {
    return std::max(std::min(x, y), v);
}

static void VS_CC dotKillInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    DotKillData *d = (DotKillData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
}

template<int sign>
static void convHoriz(const uint8_t *src, int srcStride, int16_t *dst, int width, int height) {
    while (--height) {
        dst[0] = 0;
        dst[1] = 0;

        for (int x = 2; x < width - 3; x++) {
            int16_t temp = -(src[x - 2] + src[x - 1]) + 2 * (src[x] + src[x + 1]) - (src[x + 2] + src[x + 3]);
            temp += sign * (-(src[x - 2 + srcStride] + src[x - 1 + srcStride]) + 2 * (src[x + srcStride] + src[x + 1 + srcStride]) - (src[x + 2 + srcStride] + src[x + 3 + srcStride]));
            dst[x] = temp;
        }

        dst[width - 3] = 0;
        dst[width - 2] = 0;
        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }
    memset(dst, 0, sizeof(int16_t) * width);
}

template<int sign>
static void convVert(const uint8_t *src, int srcStride, int16_t *dst, int width, int height) {
    height -= 5;

    memset(dst, 0, sizeof(int16_t) * width * 2);
    src += 2 * srcStride;
    dst += 2 * width;

    while (height--) {
        for (int x = 0; x < width - 1; x++) {
            dst[x] = -(src[x - 2 * srcStride] + src[x - 2 * srcStride + 1] + sign * (src[x - 1 * srcStride] + src[x - 1 * srcStride + 1]))
                + 2 * (src[x] + src[x + 1] + sign * (src[x + 1 * srcStride] + src[x + 1 * srcStride + 1]))
                - (src[x + 2 * srcStride] + src[x + 2 * srcStride + 1] + sign * (src[x + 3 * srcStride] + src[x + 3 * srcStride + 1]));
        }

        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }

    memset(dst, 0, sizeof(int16_t) * width * 3);
}

template<int sign>
static void applyMask(const int16_t *maskPtr, uint8_t *dst, int dstStride, int width, int height, uint8_t *ppMask) {
    maskPtr += width;
    ppMask += width;
    dst += dstStride;

    int16_t sortArray[8];

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 2; x++) {
            sortArray[0] = maskPtr[x - width - 1];
            sortArray[1] = maskPtr[x - width];
            sortArray[2] = maskPtr[x - width + 1];
            sortArray[3] = maskPtr[x - 1];
            sortArray[4] = maskPtr[x + 1];
            sortArray[5] = maskPtr[x + width - 1];
            sortArray[6] = maskPtr[x + width];
            sortArray[7] = maskPtr[x + width + 1];

            std::sort(sortArray, sortArray + 8);

            int16_t upper = sortArray[7];
            int16_t lower = sortArray[0];

            int16_t t = maskPtr[x] - clamp(maskPtr[x], lower, upper);

            if (t >= 0)
                t = (t + 4) / 8;
            else
                t = (t - 4) / 8;

            if (std::abs(t) > 1) {
                ppMask[x] = 255;

                dst[x] = (uint8_t)clamp(dst[x] - t, 0, 255);
                dst[x + 1] = (uint8_t)clamp(dst[x + 1] - t, 0, 255);
                dst[x + dstStride] = (uint8_t)clamp(dst[x + dstStride] - sign * t, 0, 255);
                dst[x + 1 + dstStride] = (uint8_t)clamp(dst[x + 1 + dstStride] - sign * t, 0, 255);
            }
        }

        maskPtr += width;
        ppMask += width;
        dst += dstStride;
    }
}

static const VSFrameRef *VS_CC dotKillGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillData *d = (DotKillData *) * instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrameRef *inframe = vsapi->getFrameFilter(n, d->node, frameCtx);
        VSFrameRef *outframe = vsapi->copyFrame(inframe, core);

        int err = 0;
        int64_t match = vsapi->propGetInt(vsapi->getFramePropsRO(inframe), "VFMMatch", 0, &err);
        vsapi->freeFrame(inframe);
        bool isC = (d->ignoreMatch || err || match == 1);
     
        int width = d->vi->width;
        int height = d->vi->height;
        int dstStride = vsapi->getStride(outframe, 0);
        uint8_t *dstPtr = vsapi->getWritePtr(outframe, 0);

        int16_t *tempMask = new int16_t[width * height];
        uint8_t *ppMask = new uint8_t[width * height]();

        for (int i = 0; i < d->napply; i++) {
            if (isC) {
                convVert<1>(dstPtr, dstStride, tempMask, width, height);
                applyMask<1>(tempMask, dstPtr, dstStride, width, height, ppMask);

                convHoriz<1>(dstPtr, dstStride, tempMask, width, height);
                applyMask<1>(tempMask, dstPtr, dstStride, width, height, ppMask);

            } else {
                convVert<-1>(dstPtr, dstStride, tempMask, width, height);
                applyMask<-1>(tempMask, dstPtr, dstStride, width, height, ppMask);

                convHoriz<-1>(dstPtr, dstStride, tempMask, width, height);
                applyMask<-1>(tempMask, dstPtr, dstStride, width, height, ppMask);
            }
        }
      
        delete[] tempMask;
        delete[] ppMask;

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    DotKillData *d = (DotKillData *)instanceData;
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillData> d(new DotKillData());

    int err;
    d->node = vsapi->propGetNode(in, "clip", 0, 0);
    d->vi = vsapi->getVideoInfo(d->node);
    d->napply = int64ToIntS(vsapi->propGetInt(in, "napply", 0, &err));
    if (d->napply < 1)
        d->napply = 1;
    d->ignoreMatch = !!vsapi->propGetInt(in, "ignorematch", 0, &err);

    vsapi->createFilter(in, out, "DotKill", dotKillInit, dotKillGetFrame, dotKillFree, fmParallel, 0, d.release(), core);
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("com.vapoursynth.dotkill", "dotkill", "VapourSynth DotKill", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("DotKill", "clip:clip;napply:int:opt;ignorematch:int:opt;", dotKillCreate, 0, plugin);
}
