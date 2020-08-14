#include "VapourSynth.h"
#include "VSHelper.h"
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <vector>

typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;
    int iterations;
    bool usematch;
} DotKillSData;


static void VS_CC dotKillSInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    DotKillSData *d = (DotKillSData *)*instanceData;
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

            int16_t t = maskPtr[x] - std::clamp(maskPtr[x], lower, upper);
            
            if (t >= 0)
                t = (t + 4) / 8;
            else
                t = (t - 4) / 8;

            if (std::abs(t) > 1) {
                ppMask[x] = 255;

                dst[x] = (uint8_t)std::clamp(dst[x] - t, 0, 255);
                dst[x + 1] = (uint8_t)std::clamp(dst[x + 1] - t, 0, 255);
                dst[x + dstStride] = (uint8_t)std::clamp(dst[x + dstStride] - sign * t, 0, 255);
                dst[x + 1 + dstStride] = (uint8_t)std::clamp(dst[x + 1 + dstStride] - sign * t, 0, 255);
            }
        }

        maskPtr += width;
        ppMask += width;
        dst += dstStride;
    }
}

static const VSFrameRef *VS_CC dotKillSGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillSData *d = (DotKillSData *)*instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrameRef *inframe = vsapi->getFrameFilter(n, d->node, frameCtx);
        VSFrameRef *outframe = vsapi->copyFrame(inframe, core);

        int err = 0;
        int64_t match = vsapi->propGetInt(vsapi->getFramePropsRO(inframe), "VFMMatch", 0, &err);
        vsapi->freeFrame(inframe);
        bool isC = (!d->usematch || err || match == 1);

        int width = d->vi->width;
        int height = d->vi->height;
        int dstStride = vsapi->getStride(outframe, 0);
        uint8_t *dstPtr = vsapi->getWritePtr(outframe, 0);

        int16_t *tempMask = new int16_t[width * height];
        uint8_t *ppMask = new uint8_t[width * height]();

        for (int i = 0; i < d->iterations; i++) {
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

static void VS_CC dotKillSFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    DotKillSData *d = (DotKillSData *)instanceData;
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillSCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillSData> d(new DotKillSData());

    int err;
    d->node = vsapi->propGetNode(in, "clip", 0, 0);
    d->vi = vsapi->getVideoInfo(d->node);
    d->iterations = int64ToIntS(vsapi->propGetInt(in, "iterations", 0, &err));
    if (d->iterations < 1)
        d->iterations = 1;
    d->usematch = !!vsapi->propGetInt(in, "usematch", 0, &err);

    vsapi->createFilter(in, out, "DotKillS", dotKillSInit, dotKillSGetFrame, dotKillSFree, fmParallel, 0, d.release(), core);
}


///////////////////////

constexpr int blockx = 16;
constexpr int blocky = 8;

// fixme, have a scene change ratio that's configurable
// if more than 3/4 of the blocks change then don't do temporal filtering at all

static void calcDiffMetric(const VSFrameRef *f1, const VSFrameRef *f2, int64_t *bdiffs, int nxblocks, int nyblocks,int field, const VSAPI *vsapi) {
    for (int plane = 0; plane < 3; plane++) {
        ptrdiff_t stride = vsapi->getStride(f1, plane);
        const uint8_t *f1p = vsapi->getReadPtr(f1, plane);
        const uint8_t *f2p = vsapi->getReadPtr(f2, plane);
        const VSFormat *fi = vsapi->getFrameFormat(f1);

        if (field) {
            f1p += stride;
            f2p += stride;
        }

        int width = vsapi->getFrameWidth(f1, plane);
        int height = vsapi->getFrameHeight(f1, plane);
        int hblockx = blockx / 2;
        int hblocky = blocky / 2;
        // adjust for subsampling
        if (plane > 0) {
            hblockx /= 1 << fi->subSamplingW;
            hblocky /= 1 << fi->subSamplingH;
        }

        for (int y = 0; y < height / 2; y++) {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx) {
                int acc = 0;
                int m = VSMIN(width, x + hblockx);
                for (int xl = x; xl < m; xl++) {
                    int tmp = f1p[xl] - f2p[xl];
                    acc += tmp * tmp;
                }
                bdiffs[ydest * nxblocks + xdest] += acc;
                xdest++;
            }

            f1p += stride * 2;
            f2p += stride * 2;
        }
    }
}

static int64_t getMaxDiff(int i, int j, const int64_t *bdiffs1, int nxblocks, int nyblocks) {
    int64_t tmp1 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j + 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j + 1];
    int64_t tmp2 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j - 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j - 1];
    int64_t tmp3 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks - j + 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j + 1];
    int64_t tmp4 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks - j - 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j - 1];
    return std::max({ tmp1, tmp2, tmp3, tmp4 });
}

static void diffMetricToMask(uint8_t *mask, const int64_t *bdiffs1, const int64_t *bdiffs2, int nxblocks, int nyblocks, int dupthresh, const VSAPI *vsapi) {
    // fixme, must fill in edge blocks
    int totdiff1 = 0;
    int totdiff2 = 0;

    for (int i = 1; i < nyblocks - 1; i++) {
        for (int j = 1; j < nxblocks - 1; j++) {
            int64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks, nyblocks);
            int64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks, nyblocks);

            if (diff1 >= dupthresh)
                totdiff1++;
            if (diff2 >= dupthresh)
                totdiff2++;
        }
    }

    // 1/3 of the blocks are different = no temporal processing since it's most likely a pan that covers most of the screen
    bool skip1 = (totdiff1 * 3 > nxblocks * nyblocks);
    bool skip2 = (totdiff2 * 3 > nxblocks * nyblocks);

    for (int i = 1; i < nyblocks - 1; i++) {
        for (int j = 1; j < nxblocks - 1; j++) {
            int64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks, nyblocks);
            int64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks, nyblocks);

            if (!skip1 && diff1 <= diff2 && diff1 < dupthresh)
                mask[nxblocks * i + j] = 1;
            else if (!skip2 && diff2 < diff1 && diff2 < dupthresh)
                mask[nxblocks * i + j] = 2;
            else
                mask[nxblocks * i + j] = 0;
        }

        // extend mask left and right
        mask[nxblocks * i] = mask[nxblocks * i + 1];
        mask[nxblocks * i + (nxblocks - 1)] = mask[nxblocks * i + (nxblocks - 2)];
    }

    // extend mask to top and bottom
    memcpy(mask, mask + nxblocks, nxblocks);
    memcpy(mask + nxblocks * (nyblocks - 1), mask + nxblocks * (nyblocks - 2), nxblocks);
}

static int64_t calcMetric(const VSFrameRef *f1, const VSFrameRef *f2, uint8_t *mask, int nxblocks, int nyblocks, int64_t &maxdiff, int field, int dupthresh, const VSAPI *vsapi) {
    std::vector<int64_t> bdiffs(nxblocks * nyblocks);
    for (int plane = 0; plane < 3; plane++) {
        ptrdiff_t stride = vsapi->getStride(f1, plane);
        const uint8_t *f1p = vsapi->getReadPtr(f1, plane);
        const uint8_t *f2p = vsapi->getReadPtr(f2, plane);
        const VSFormat *fi = vsapi->getFrameFormat(f1);

        if (field) {
            f1p += stride;
            f2p += stride;
        }

        int width = vsapi->getFrameWidth(f1, plane);
        int height = vsapi->getFrameHeight(f1, plane);
        int hblockx = blockx / 2;
        int hblocky = blocky / 2;
        // adjust for subsampling
        if (plane > 0) {
            hblockx /= 1 << fi->subSamplingW;
            hblocky /= 1 << fi->subSamplingH;
        }

        for (int y = 0; y < height / 2; y++) {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx) {
                int acc = 0;
                int m = VSMIN(width, x + hblockx);
                for (int xl = x; xl < m; xl++) {
                    int tmp = f1p[xl] - f2p[xl];
                    acc += tmp * tmp;
                }
                bdiffs[ydest * nxblocks + xdest] += acc;
                xdest++;
            }

            f1p += stride * 2;
            f2p += stride * 2;
        }
    }

    // fixme, edge blocks are unset
    for (int i = 1; i < nyblocks - 1; i++) {
        for (int j = 1; j < nxblocks - 1; j++) {
            int64_t tmp1 = bdiffs[i * nxblocks + j] + bdiffs[i * nxblocks + j + 1] + bdiffs[(i + 1) * nxblocks + j] + bdiffs[(i + 1) * nxblocks + j + 1];
            int64_t tmp2 = bdiffs[i * nxblocks + j] + bdiffs[i * nxblocks + j - 1] + bdiffs[(i + 1) * nxblocks + j] + bdiffs[(i + 1) * nxblocks + j - 1];
            int64_t tmp3 = bdiffs[i * nxblocks + j] + bdiffs[i * nxblocks - j + 1] + bdiffs[(i - 1) * nxblocks + j] + bdiffs[(i - 1) * nxblocks + j + 1];
            int64_t tmp4 = bdiffs[i * nxblocks + j] + bdiffs[i * nxblocks - j - 1] + bdiffs[(i - 1) * nxblocks + j] + bdiffs[(i - 1) * nxblocks + j - 1];
            tmp1 = std::max({ tmp1, tmp2, tmp3, tmp4 });
            mask[nxblocks * i + j] = tmp1 < dupthresh ? 255 : 0;
        }
    }

    int64_t maxdiff1 = -1;
    for (int i = 0; i < nyblocks - 1; i++) {
        for (int j = 0; j < nxblocks - 1; j++) {
            int64_t tmp = bdiffs[i * nxblocks + j] + bdiffs[i * nxblocks + j + 1] + bdiffs[(i + 1) * nxblocks + j] + bdiffs[(i + 1) * nxblocks + j + 1];
            if (tmp > maxdiff1)
                maxdiff1 = tmp;
        }
    }
    maxdiff = maxdiff1;

    int64_t totdiff1 = 0;
    for (int i = 0; i < nxblocks * nyblocks; i++) {
        // saturate on overflow?
        assert(totdiff1 + bdiffs[i] >= totdiff1);
        totdiff1 += bdiffs[i];
    }

    return totdiff1;
}

typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;
    int order;
    int offset;
} DotKillZData;


static void VS_CC dotKillZInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    DotKillZData *d = (DotKillZData *)*instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
}


static void applyStaticMask2(VSFrameRef *dst, const VSFrameRef *f0, const VSFrameRef *f1, const VSFrameRef *f2, uint8_t *mask, int nxblocks, int nyblocks, int field, const VSAPI *vsapi) {
    for (int plane = 0; plane < 3; plane++) {
        ptrdiff_t stride = vsapi->getStride(f1, plane);
        const uint8_t *f0p = vsapi->getReadPtr(f0, plane);
        const uint8_t *f1p = vsapi->getReadPtr(f1, plane);
        const uint8_t *f2p = vsapi->getReadPtr(f2, plane);
        uint8_t *dstp = vsapi->getWritePtr(dst, plane);
        const VSFormat *fi = vsapi->getFrameFormat(f1);

        if (field) {
            f0p += stride;
            f1p += stride;
            f2p += stride;
            dstp += stride;
        }

        int width = vsapi->getFrameWidth(f1, plane);
        int height = vsapi->getFrameHeight(f1, plane);
        int hblockx = blockx / 2;
        int hblocky = blocky / 2;
        // adjust for subsampling
        if (plane > 0) {
            hblockx /= 1 << fi->subSamplingW;
            hblocky /= 1 << fi->subSamplingH;
        }

        for (int y = 0; y < height / 2; y++) {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx) {
                int acc = 0;
                int m = VSMIN(width, x + hblockx);

                for (int xl = x; xl < m; xl++) {
                    if (mask[ydest * nxblocks + xdest] == 1)
                        dstp[xl] = ((f1p[xl] + f0p[xl] + 1) / 2);
                    else if (mask[ydest * nxblocks + xdest] == 2)
                        dstp[xl] = ((f2p[xl] + f0p[xl] + 1) / 2);
                    else
                        //dstp[xl] = 255;
                        dstp[xl] = dstp[xl];
                }
                xdest++;
            }

            f0p += stride * 2;
            f1p += stride * 2;
            f2p += stride * 2;
            dstp += stride * 2;
        }
    }
}

static void applyFieldBlend(const VSFrameRef *srcc, const VSFrameRef *srcn, VSFrameRef *outframe, int order, VSCore *core, const VSAPI *vsapi) {
    const VSFormat *fi = vsapi->getFrameFormat(srcc);
    for (int plane = 0; plane < fi->numPlanes; plane++) {
        int width = vsapi->getFrameWidth(srcc, plane);
        int height = vsapi->getFrameHeight(srcc, plane);

        int stride = vsapi->getStride(outframe, plane);
        uint8_t *dstp = vsapi->getWritePtr(outframe, plane);
        const uint8_t *srccp = vsapi->getReadPtr(srcc, plane);
        const uint8_t *srcnp = vsapi->getReadPtr(srcn, plane);

        if (order) {
            srccp += stride;
            srcnp += stride;
            dstp += stride;
        }

        for (int h = order; h < height; h += 2) {
            for (int w = 0; w < width; w++)
                dstp[w] = (srccp[w] + srcnp[w] + 1) / 2;

            srccp += 2 * stride;
            srcnp += 2 * stride;
            dstp += 2 * stride;
        }
    }
}

static void applyDotcrawInverse(const VSFrameRef *srcc, const VSFrameRef *srcn, VSFrameRef *outframe, int order, VSCore *core, const VSAPI *vsapi) {
    const VSFormat *fi = vsapi->getFrameFormat(srcc);
    for (int plane = 0; plane < fi->numPlanes; plane++) {
        int width = vsapi->getFrameWidth(srcc, plane);
        int height = vsapi->getFrameHeight(srcc, plane);

        int stride = vsapi->getStride(outframe, plane);
        uint8_t *dstp = vsapi->getWritePtr(outframe, plane);
        const uint8_t *srccp = vsapi->getReadPtr(srcc, plane);
        const uint8_t *srcnp = vsapi->getReadPtr(srcn, plane);

        if (order) {
            srccp += stride;
            srcnp += stride;
            dstp += stride;
        }

        for (int h = order; h < height; h += 2) {
            for (int w = 0; w < width; w++) {
                dstp[w] = (srccp[w] + srcnp[w] + 1) / 2;

                if (h > 1) {
                    uint8_t l0val = dstp[w - 2 * stride];
                    uint8_t l2val = dstp[w];
                    int l0diff = dstp[w - 2 * stride] - srccp[w - 2 * stride];
                    int l2diff = dstp[w] - srccp[w];
                    if (plane == 0)
                        dstp[w - stride] = std::clamp(srccp[w - stride] + (order ? l0diff : l2diff), 16, 235);
                    else
                        dstp[w - stride] = (l0val + l2val + 1) / 2; // simply use some kind of interpolation and discard one field?
                }
            }

            srccp += 2 * stride;
            srcnp += 2 * stride;
            dstp += 2 * stride;
        }
    }
}

static const VSFrameRef *VS_CC dotKillZGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillZData *d = (DotKillZData *)*instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrameRef *srcp = vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        const VSFrameRef *srcc = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef *srcn = vsapi->getFrameFilter(n + 1, d->node, frameCtx);

        /*
            FIELD OFFSETS
            -1  0  1  2  3
            A1 B1 B1 C1 D1
            A2 B2 C2 D2 D2
        */

        VSFrameRef *outframe = vsapi->copyFrame(srcc, core);
        if ((n + d->offset) % 5 == 0) {
            // current and next field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
            applyDotcrawInverse(srcc, srcn, outframe, d->order, core, vsapi);
        } else if ((n + d->offset) % 5 == 1) {
            // current and previous field are duplicates so blend them together
            applyFieldBlend(srcc, srcp, outframe, d->order, core, vsapi);
        } else if ((n + d->offset) % 5 == 2) {
            // current and next complement field are duplicates so blend them together
            applyFieldBlend(srcc, srcn, outframe, !d->order, core, vsapi);
        } else if ((n + d->offset) % 5 == 3) {
            // current and previous field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
            applyDotcrawInverse(srcc, srcp, outframe, !d->order, core, vsapi);
        }

        vsapi->freeFrame(srcp);
        vsapi->freeFrame(srcc);
        vsapi->freeFrame(srcn);

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillZFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    DotKillZData *d = (DotKillZData *)instanceData;
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillZCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillZData> d(new DotKillZData());

    int err;
    d->node = vsapi->propGetNode(in, "clip", 0, 0);
    d->vi = vsapi->getVideoInfo(d->node);
    d->offset = int64ToIntS(vsapi->propGetInt(in, "offset", 0, &err));
    d->order = !!vsapi->propGetInt(in, "order", 0, &err);

    vsapi->createFilter(in, out, "DotKillZ", dotKillZInit, dotKillZGetFrame, dotKillZFree, fmParallelRequests, 0, d.release(), core);
}

//////////////
// DotKillT

typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;
    int order;
    int offset;
    int dupthresh;
} DotKillTData;


static const VSFrameRef *VS_CC dotKillTGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillTData *d = (DotKillTData *)*instanceData;

    int nxblocks = (d->vi->width + blockx / 2 - 1) / (blockx / 2);
    int nyblocks = (d->vi->height + blocky / 2 - 1) / (blocky / 2);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(std::max(n - 2, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 2, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrameRef *srcpp = vsapi->getFrameFilter(std::max(n - 2, 0), d->node, frameCtx);
        const VSFrameRef *srcp = vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        const VSFrameRef *srcc = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef *srcn = vsapi->getFrameFilter(n + 1, d->node, frameCtx);
        const VSFrameRef *srcnn = vsapi->getFrameFilter(n + 2, d->node, frameCtx);

        // first two fields are duplicates, meaning that offset 0, 1 and 5, 6 are used in the various calculations
        // fields 2, 3 and 4 needs to have determined how many blocks are consecutively static
        // note that comparisons are only run on a single since we can in most comparisons can eliminate
        // the dotcrawl from the equation

        // 0-1 2 3 4 5-6

        // the complementary field can likewise be used for movement detection

        /*
            FIELD OFFSETS
            -1  0  1  2  3  4
            A1 B1 B1 C1 D1 E1
            A2 B2 C2 D2 D2 E2
        */

        VSFrameRef *outframe = vsapi->copyFrame(srcc, core);

        vsapi->propSetInt(vsapi->getFramePropsRW(outframe), "offset", (n + d->offset) % 5, paReplace);

        int64_t maxdiff = 0;

        if ((n + d->offset) % 5 == 0) {
            // 1
            applyDotcrawInverse(srcc, srcn, outframe, d->order, core, vsapi);

            // 2
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric(srcp, srcn, maskprev1.data(), nxblocks, nyblocks, d->order, vsapi);
            calcDiffMetric(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

            applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !d->order, vsapi);
        } else if ((n + d->offset) % 5 == 1) {
            // 1
            applyFieldBlend(srcc, srcp, outframe, d->order, core, vsapi);

            // 2
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric(srcp, srcn, maskprev1.data(), nxblocks, nyblocks, d->order, vsapi);
            calcDiffMetric(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, !d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

            applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !d->order, vsapi);
        } else if ((n + d->offset) % 5 == 2) {
            // 1
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, d->order, vsapi);
            calcDiffMetric(srcp, srcn, masknext1.data(), nxblocks, nyblocks, !d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

            applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, d->order, vsapi);

            // 2
            applyFieldBlend(srcc, srcn, outframe, !d->order, core, vsapi);
        } else if ((n + d->offset) % 5 == 3) {
            // 2
            applyDotcrawInverse(srcc, srcp, outframe, !d->order, core, vsapi);

            // 1
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, !d->order, vsapi);
            calcDiffMetric(srcp, srcn, masknext1.data(), nxblocks, nyblocks, !d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

            applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, d->order, vsapi);
        } else if ((n + d->offset) % 5 == 4) {
            // 1
            {
                std::vector<int64_t> maskprev1(nxblocks * nyblocks);
                std::vector<int64_t> masknext1(nxblocks * nyblocks);
                std::vector<uint8_t> mask(nxblocks * nyblocks);
                calcDiffMetric(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, !d->order, vsapi);
                calcDiffMetric(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, d->order, vsapi);

                diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

                applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, d->order, vsapi);
            }

            // 2
            {
                std::vector<int64_t> maskprev1(nxblocks * nyblocks);
                std::vector<int64_t> masknext1(nxblocks * nyblocks);
                std::vector<uint8_t> mask(nxblocks * nyblocks);
                calcDiffMetric(srcpp, srcc, maskprev1.data(), nxblocks, nyblocks, !d->order, vsapi);
                calcDiffMetric(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, d->order, vsapi);

                diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, vsapi);

                applyStaticMask2(outframe, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !d->order, vsapi);
            }
        }

        vsapi->freeFrame(srcpp);
        vsapi->freeFrame(srcp);
        vsapi->freeFrame(srcc);
        vsapi->freeFrame(srcn);
        vsapi->freeFrame(srcnn);

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillTFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    DotKillTData *d = (DotKillTData *)instanceData;
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillTInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    DotKillTData *d = (DotKillTData *)*instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
}


static void VS_CC dotKillTCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillTData> d(new DotKillTData());

    int err;
    d->node = vsapi->propGetNode(in, "clip", 0, 0);
    d->vi = vsapi->getVideoInfo(d->node);
    d->offset = int64ToIntS(vsapi->propGetInt(in, "offset", 0, &err));
    d->order = !!vsapi->propGetInt(in, "order", 0, &err);
    d->dupthresh = int64ToIntS(vsapi->propGetInt(in, "dupthresh", 0, &err));
    if (err || d->dupthresh < 0)
        d->dupthresh = 64;
    d->dupthresh *= d->dupthresh;

    vsapi->createFilter(in, out, "DotKillT", dotKillTInit, dotKillTGetFrame, dotKillTFree, fmParallelRequests, 0, d.release(), core);
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("com.vapoursynth.dotkill", "dotkill", "VapourSynth DotKill", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("DotKillS", "clip:clip;iterations:int:opt;usematch:int:opt;", dotKillSCreate, 0, plugin);
    registerFunc("DotKillZ", "clip:clip;order:int:opt;offset:int:opt;", dotKillZCreate, 0, plugin);
    registerFunc("DotKillT", "clip:clip;order:int:opt;offset:int:opt;dupthresh:int:opt;", dotKillTCreate, 0, plugin);
}
