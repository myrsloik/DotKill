#include "VapourSynth4.h"
#include "VSHelper4.h"
#include <algorithm>
#include <memory>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <vector>

// The convolutions have a gain of 8 and a sign bit, so the intermediate needs 4 bits
// more than the sample. That still fits int16_t for 8 bit input, which keeps the
// temporary buffer and its memory traffic half the size of the 16 bit one.
template<typename T>
using mask_t = std::conditional_t<sizeof(T) == 1, int16_t, int32_t>;

// Squared sample differences accumulated over one block row. A single squared 16 bit
// difference is already 65535*65535, which overflows int32_t on its own, while the 8
// bit worst case of 8 taps times 255*255 stays comfortably inside it.
template<typename T>
using diff_t = std::conditional_t<sizeof(T) == 1, int, int64_t>;

// VapourSynth reports strides in bytes while the processing indexes whole samples
template<typename T>
static ptrdiff_t elementStride(const VSFrame *f, int plane, const VSAPI *vsapi) {
    return vsapi->getStride(f, plane) / static_cast<ptrdiff_t>(sizeof(T));
}

// Every integer depth from 8 to 16 is handled by the uint8_t and uint16_t
// instantiations. VapourSynth also permits 17 to 32 bit integer, where
// bytesPerSample becomes 4, so the upper bound here is what keeps the
// bytesPerSample based dispatch in the getframe functions honest.
static bool isSupportedFormat(const VSVideoInfo *vi, bool allowGray) {
    if (!vsh::isConstantVideoFormat(vi))
        return false;

    const VSVideoFormat &fi = vi->format;
    if (fi.sampleType != stInteger || fi.bitsPerSample < 8 || fi.bitsPerSample > 16)
        return false;

    // Pinning the subsampling to 420 is load bearing beyond matching the algorithm:
    // calcDiffMetric and applyTemporalMask derive the chroma block step by dividing
    // blockx/2 and blocky/2 by the subsampling, and a wider ratio would round that
    // step down to zero, turning the x loop into an endless one.
    if (fi.colorFamily == cfYUV)
        return fi.subSamplingW == 1 && fi.subSamplingH == 1;

    return allowGray && fi.colorFamily == cfGray;
}

////////////////////////////////////////
// DotKillS

typedef struct {
    VSNode *node;
    const VSVideoInfo *vi;
    int iterations;
} DotKillSData;

template<typename T>
static void convHoriz(const T *src, ptrdiff_t srcStride, mask_t<T> *dst, int width, int height) {
    while (--height) {
        dst[0] = 0;
        dst[1] = 0;

        for (int x = 2; x < width - 3; x++) {
            mask_t<T> temp = -(src[x - 2] + src[x - 1]) + 2 * (src[x] + src[x + 1]) - (src[x + 2] + src[x + 3]);
            temp += -(src[x - 2 + srcStride] + src[x - 1 + srcStride]) + 2 * (src[x + srcStride] + src[x + 1 + srcStride]) - (src[x + 2 + srcStride] + src[x + 3 + srcStride]);
            dst[x] = temp;
        }

        dst[width - 3] = 0;
        dst[width - 2] = 0;
        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }
    std::memset(dst, 0, sizeof(*dst) * width);
}

template<typename T>
static void convVert(const T *src, ptrdiff_t srcStride, mask_t<T> *dst, int width, int height) {
    height -= 5;

    std::memset(dst, 0, sizeof(*dst) * width * 2);
    src += 2 * srcStride;
    dst += 2 * width;

    while (height--) {
        for (int x = 0; x < width - 1; x++) {
            dst[x] = -(src[x - 2 * srcStride] + src[x - 2 * srcStride + 1] + (src[x - 1 * srcStride] + src[x - 1 * srcStride + 1]))
                + 2 * (src[x] + src[x + 1] + (src[x + 1 * srcStride] + src[x + 1 * srcStride + 1]))
                - (src[x + 2 * srcStride] + src[x + 2 * srcStride + 1] + (src[x + 3 * srcStride] + src[x + 3 * srcStride + 1]));
        }

        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }

    std::memset(dst, 0, sizeof(*dst) * width * 3);
}

template<typename T>
static void applyMask(const mask_t<T> *maskPtr, T *dst, ptrdiff_t dstStride, int width, int height, int shift) {
    const int rangeLow = 16 << shift;
    const int rangeHigh = 235 << shift;
    // the minimum correction worth applying is expressed in 8 bit units so it has
    // to track the sample range, otherwise deeper input filters far more eagerly
    const mask_t<T> threshold = static_cast<mask_t<T>>(2) << shift;

    maskPtr += width;
    dst += dstStride;

    mask_t<T> sortArray[8];

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

            mask_t<T> upper = sortArray[7];
            mask_t<T> lower = sortArray[0];

            mask_t<T> t = maskPtr[x] - std::clamp(maskPtr[x], lower, upper);

            if (t >= 0)
                t = (t + 4) / 8;
            else
                t = (t - 4) / 8;

            if (std::abs(t) >= threshold) {
                dst[x] = static_cast<T>(std::clamp<int>(dst[x] - t, rangeLow, rangeHigh));
                dst[x + 1] = static_cast<T>(std::clamp<int>(dst[x + 1] - t, rangeLow, rangeHigh));
                dst[x + dstStride] = static_cast<T>(std::clamp<int>(dst[x + dstStride] - t, rangeLow, rangeHigh));
                dst[x + 1 + dstStride] = static_cast<T>(std::clamp<int>(dst[x + 1 + dstStride] - t, rangeLow, rangeHigh));
            }
        }

        maskPtr += width;
        dst += dstStride;
    }
}

template<typename T>
static void dotKillSProcess(VSFrame *outframe, int width, int height, int shift, int iterations, const VSAPI *vsapi) {
    ptrdiff_t dstStride = elementStride<T>(outframe, 0, vsapi);
    T *dstPtr = reinterpret_cast<T *>(vsapi->getWritePtr(outframe, 0));

    std::vector<mask_t<T>> tempMask(static_cast<size_t>(width) * height);

    for (int i = 0; i < iterations; i++) {
        convVert<T>(dstPtr, dstStride, tempMask.data(), width, height);
        applyMask<T>(tempMask.data(), dstPtr, dstStride, width, height, shift);

        convHoriz<T>(dstPtr, dstStride, tempMask.data(), width, height);
        applyMask<T>(tempMask.data(), dstPtr, dstStride, width, height, shift);
    }
}

static const VSFrame *VS_CC dotKillSGetFrame(int n, int activationReason, void *instanceData, [[maybe_unused]] void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillSData *d = reinterpret_cast<DotKillSData*>(instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *inframe = vsapi->getFrameFilter(n, d->node, frameCtx);
        VSFrame *outframe = vsapi->copyFrame(inframe, core);
        vsapi->freeFrame(inframe);

        int shift = d->vi->format.bitsPerSample - 8;

        if (d->vi->format.bytesPerSample == 1)
            dotKillSProcess<uint8_t>(outframe, d->vi->width, d->vi->height, shift, d->iterations, vsapi);
        else
            dotKillSProcess<uint16_t>(outframe, d->vi->width, d->vi->height, shift, d->iterations, vsapi);

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillSFree(void *instanceData, [[maybe_unused]] VSCore *core, const VSAPI *vsapi) {
    DotKillSData *d = reinterpret_cast<DotKillSData*>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillSCreate(const VSMap *in, VSMap *out, [[maybe_unused]] void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillSData> d(new DotKillSData());

    int err;
    d->node = vsapi->mapGetNode(in, "clip", 0, nullptr);
    d->vi = vsapi->getVideoInfo(d->node);
    d->iterations = std::clamp(vsapi->mapGetIntSaturated(in, "iterations", 0, &err), 1, 10);
    if (!isSupportedFormat(d->vi, true)) {
        vsapi->mapSetError(out, "DotKillS: only constant dimension 8-16 bit integer YUV420 and GRAY supported");
        vsapi->freeNode(d->node);
        return;
    }

    // convHoriz taps x-2..x+3 and convVert taps 3 rows in both directions, so
    // anything smaller than the kernel support underflows the temporary buffer
    if (d->vi->width < 6 || d->vi->height < 6) {
        vsapi->mapSetError(out, "DotKillS: clip dimensions must be at least 6x6");
        vsapi->freeNode(d->node);
        return;
    }

    VSFilterDependency deps[] = { {d->node, rpStrictSpatial} };
    vsapi->createVideoFilter(out, "DotKillS", d->vi, dotKillSGetFrame, dotKillSFree, fmParallel, deps, 1, d.get(), core);
    d.release();
}

/////////////////////////////////////////////////////////////////////
// DotKillZ


template<typename T>
static void applyFieldBlend(const VSFrame *srcc, const VSFrame *srcn, VSFrame *outframe, int order, const VSAPI *vsapi) {
    const VSVideoFormat *fi = vsapi->getVideoFrameFormat(srcc);
    for (int plane = 0; plane < fi->numPlanes; plane++) {
        int width = vsapi->getFrameWidth(srcc, plane);
        int height = vsapi->getFrameHeight(srcc, plane);

        ptrdiff_t dstStride = elementStride<T>(outframe, plane, vsapi);
        ptrdiff_t srccStride = elementStride<T>(srcc, plane, vsapi);
        ptrdiff_t srcnStride = elementStride<T>(srcn, plane, vsapi);
        T *dstp = reinterpret_cast<T *>(vsapi->getWritePtr(outframe, plane));
        const T *srccp = reinterpret_cast<const T *>(vsapi->getReadPtr(srcc, plane));
        const T *srcnp = reinterpret_cast<const T *>(vsapi->getReadPtr(srcn, plane));

        if (order) {
            srccp += srccStride;
            srcnp += srcnStride;
            dstp += dstStride;
        }

        for (int h = order; h < height; h += 2) {
            for (int w = 0; w < width; w++)
                dstp[w] = static_cast<T>((srccp[w] + srcnp[w] + 1) / 2);

            srccp += 2 * srccStride;
            srcnp += 2 * srcnStride;
            dstp += 2 * dstStride;
        }
    }
}

template<typename T>
static void applyDotcrawInverse(const VSFrame *srcc, const VSFrame *srcn, VSFrame *outframe, int order, int shift, const VSAPI *vsapi) {
    const int rangeLow = 16 << shift;
    const int rangeHigh = 235 << shift;

    const VSVideoFormat *fi = vsapi->getVideoFrameFormat(srcc);
    for (int plane = 0; plane < fi->numPlanes; plane++) {
        int width = vsapi->getFrameWidth(srcc, plane);
        int height = vsapi->getFrameHeight(srcc, plane);

        ptrdiff_t dstStride = elementStride<T>(outframe, plane, vsapi);
        ptrdiff_t srccStride = elementStride<T>(srcc, plane, vsapi);
        ptrdiff_t srcnStride = elementStride<T>(srcn, plane, vsapi);
        T *dstp = reinterpret_cast<T *>(vsapi->getWritePtr(outframe, plane));
        const T *srccp = reinterpret_cast<const T *>(vsapi->getReadPtr(srcc, plane));
        const T *srcnp = reinterpret_cast<const T *>(vsapi->getReadPtr(srcn, plane));

        if (order) {
            srccp += srccStride;
            srcnp += srcnStride;
            dstp += dstStride;
        }

        for (int h = order; h < height; h += 2) {
            for (int w = 0; w < width; w++) {
                dstp[w] = static_cast<T>((srccp[w] + srcnp[w] + 1) / 2);

                if (h > 1) {
                    T l0val = dstp[w - 2 * dstStride];
                    T l2val = dstp[w];
                    int l0diff = dstp[w - 2 * dstStride] - srccp[w - 2 * srccStride];
                    int l2diff = dstp[w] - srccp[w];
                    if (plane == 0)
                        dstp[w - dstStride] = static_cast<T>(std::clamp<int>(srccp[w - srccStride] + (order ? l0diff : l2diff), rangeLow, rangeHigh));
                    else
                        dstp[w - dstStride] = static_cast<T>((l0val + l2val + 1) / 2); // simply use some kind of interpolation and discard one field?
                }
            }

            srccp += 2 * srccStride;
            srcnp += 2 * srcnStride;
            dstp += 2 * dstStride;
        }
    }
}

typedef struct {
    VSNode *node;
    const VSVideoInfo *vi;
    int order;
    int offset;
} DotKillZData;

template<typename T>
static void dotKillZProcess(VSFrame *outframe, const VSFrame *srcp, const VSFrame *srcc, const VSFrame *srcn, int pattern, int order, int shift, const VSAPI *vsapi) {
    if (pattern == 0) {
        // current and next field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
        applyDotcrawInverse<T>(srcc, srcn, outframe, order, shift, vsapi);
    } else if (pattern == 1) {
        // current and previous field are duplicates so blend them together
        applyFieldBlend<T>(srcc, srcp, outframe, order, vsapi);
    } else if (pattern == 2) {
        // current and next complement field are duplicates so blend them together
        applyFieldBlend<T>(srcc, srcn, outframe, !order, vsapi);
    } else if (pattern == 3) {
        // current and previous field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
        applyDotcrawInverse<T>(srcc, srcp, outframe, !order, shift, vsapi);
    }
}

static const VSFrame *VS_CC dotKillZGetFrame(int n, int activationReason, void *instanceData, [[maybe_unused]] void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillZData *d = reinterpret_cast<DotKillZData*>(instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *srcp = vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        const VSFrame *srcc = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *srcn = vsapi->getFrameFilter(n + 1, d->node, frameCtx);

        /*
            FIELD OFFSETS
            -1  0  1  2  3
            A1 B1 B1 C1 D1
            A2 B2 C2 D2 D2
        */

        VSFrame *outframe = vsapi->copyFrame(srcc, core);
        int shift = d->vi->format.bitsPerSample - 8;

        if (d->vi->format.bytesPerSample == 1)
            dotKillZProcess<uint8_t>(outframe, srcp, srcc, srcn, (n + d->offset) % 5, d->order, shift, vsapi);
        else
            dotKillZProcess<uint16_t>(outframe, srcp, srcc, srcn, (n + d->offset) % 5, d->order, shift, vsapi);

        vsapi->freeFrame(srcp);
        vsapi->freeFrame(srcc);
        vsapi->freeFrame(srcn);

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillZFree(void *instanceData, [[maybe_unused]] VSCore *core, const VSAPI *vsapi) {
    DotKillZData *d = reinterpret_cast<DotKillZData*>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillZCreate(const VSMap *in, VSMap *out, [[maybe_unused]] void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillZData> d(new DotKillZData());

    int err;
    d->node = vsapi->mapGetNode(in, "clip", 0, nullptr);
    d->vi = vsapi->getVideoInfo(d->node);
    d->offset = vsapi->mapGetIntSaturated(in, "offset", 0, &err);
    d->order = !!vsapi->mapGetInt(in, "order", 0, &err);

    if (!isSupportedFormat(d->vi, false)) {
        vsapi->mapSetError(out, "DotKillZ: only constant dimension 8-16 bit integer YUV420 supported");
        vsapi->freeNode(d->node);
        return;
    }

    if (d->offset < 0 || d->offset > 4) {
        vsapi->mapSetError(out, "DotKillZ: offset must be between 0 and 4");
        vsapi->freeNode(d->node);
        return;
    }

    VSFilterDependency deps[] = { {d->node, rpGeneral} };
    vsapi->createVideoFilter(out, "DotKillZ", d->vi, dotKillZGetFrame, dotKillZFree, fmParallel, deps, 1, d.get(), core);
    d.release();
}

/////////////////////////////////////////////////////////////////////
// DotKillT

constexpr int blockx = 16;
constexpr int blocky = 8;

template<typename T>
static void calcDiffMetric(const VSFrame *f1, const VSFrame *f2, int64_t *bdiffs, int nxblocks, int field, const VSAPI *vsapi) {
    for (int plane = 0; plane < 3; plane++) {
        ptrdiff_t f1Stride = elementStride<T>(f1, plane, vsapi);
        ptrdiff_t f2Stride = elementStride<T>(f2, plane, vsapi);
        const T *f1p = reinterpret_cast<const T *>(vsapi->getReadPtr(f1, plane));
        const T *f2p = reinterpret_cast<const T *>(vsapi->getReadPtr(f2, plane));
        const VSVideoFormat *fi = vsapi->getVideoFrameFormat(f1);

        if (field) {
            f1p += f1Stride;
            f2p += f2Stride;
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
                diff_t<T> acc = 0;
                int m = VSMIN(width, x + hblockx);
                for (int xl = x; xl < m; xl++) {
                    diff_t<T> tmp = f1p[xl] - f2p[xl];
                    acc += tmp * tmp;
                }
                bdiffs[ydest * nxblocks + xdest] += acc;
                xdest++;
            }

            f1p += f1Stride * 2;
            f2p += f2Stride * 2;
        }
    }
}

static int64_t getMaxDiff(int i, int j, const int64_t *bdiffs1, int nxblocks) {
    int64_t tmp1 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j + 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j + 1];
    int64_t tmp2 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j - 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j - 1];
    int64_t tmp3 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j + 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j + 1];
    int64_t tmp4 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j - 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j - 1];
    return std::max({ tmp1, tmp2, tmp3, tmp4 });
}

static void diffMetricToMask(uint8_t *mask, const int64_t *bdiffs1, const int64_t *bdiffs2, int nxblocks, int nyblocks, int64_t dupthresh, int tratio) {
    int totdiff1 = 0;
    int totdiff2 = 0;

    for (int i = 1; i < nyblocks - 1; i++) {
        for (int j = 1; j < nxblocks - 1; j++) {
            int64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks);
            int64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks);

            if (diff1 >= dupthresh)
                totdiff1++;
            if (diff2 >= dupthresh)
                totdiff2++;
        }
    }

    // skip temporal processing if more than 1/tratio blocks have changed.
    // widened because tratio has no upper bound while the interior block count
    // reaches five figures, so the 32 bit product would overflow
    int64_t interior = static_cast<int64_t>(nxblocks - 2) * (nyblocks - 2);
    bool skip1 = (static_cast<int64_t>(totdiff1) * tratio > interior);
    bool skip2 = (static_cast<int64_t>(totdiff2) * tratio > interior);

    for (int i = 1; i < nyblocks - 1; i++) {
        for (int j = 1; j < nxblocks - 1; j++) {
            int64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks);
            int64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks);

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
    std::memcpy(mask, mask + nxblocks, nxblocks);
    std::memcpy(mask + nxblocks * (nyblocks - 1), mask + nxblocks * (nyblocks - 2), nxblocks);
}

template<typename T>
static void applyTemporalMask(VSFrame *dst, const VSFrame *f0, const VSFrame *f1, const VSFrame *f2, const uint8_t *mask, int nxblocks, int field, bool show, int shift, const VSAPI *vsapi) {
    for (int plane = 0; plane < 3; plane++) {
        ptrdiff_t dstStride = elementStride<T>(dst, plane, vsapi);
        ptrdiff_t f0Stride = elementStride<T>(f0, plane, vsapi);
        ptrdiff_t f1Stride = elementStride<T>(f1, plane, vsapi);
        ptrdiff_t f2Stride = elementStride<T>(f2, plane, vsapi);
        const T *f0p = reinterpret_cast<const T *>(vsapi->getReadPtr(f0, plane));
        const T *f1p = reinterpret_cast<const T *>(vsapi->getReadPtr(f1, plane));
        const T *f2p = reinterpret_cast<const T *>(vsapi->getReadPtr(f2, plane));
        T *dstp = reinterpret_cast<T *>(vsapi->getWritePtr(dst, plane));
        const VSVideoFormat *fi = vsapi->getVideoFrameFormat(f1);

        if (field) {
            f0p += f0Stride;
            f1p += f1Stride;
            f2p += f2Stride;
            dstp += dstStride;
        }

        int width = vsapi->getFrameWidth(f1, plane);
        int height = vsapi->getFrameHeight(f1, plane);
        int hblockx = blockx / 2;
        int hblocky = blocky / 2;

        if (plane > 0) {
            hblockx /= 1 << fi->subSamplingW;
            hblocky /= 1 << fi->subSamplingH;
        }

        for (int y = 0; y < height / 2; y++) {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx) {
                int m = VSMIN(width, x + hblockx);

                for (int xl = x; xl < m; xl++) {
                    if (mask[ydest * nxblocks + xdest] == 1)
                        dstp[xl] = static_cast<T>((f1p[xl] + f0p[xl] + 1) / 2);
                    else if (mask[ydest * nxblocks + xdest] == 2)
                        dstp[xl] = static_cast<T>((f2p[xl] + f0p[xl] + 1) / 2);
                }
                xdest++;
            }

            f0p += f0Stride * 2;
            f1p += f1Stride * 2;
            f2p += f2Stride * 2;
            dstp += dstStride * 2;
        }
    }

    // Horrible square drawing code
    if (show) {
        // the markers are full scale rather than TV range, same as the 8 bit original
        const T markLow = 0;
        const T markHigh = static_cast<T>((1 << (8 + shift)) - 1);

        ptrdiff_t stride = elementStride<T>(dst, 0, vsapi);
        T *dstp = reinterpret_cast<T *>(vsapi->getWritePtr(dst, 0));

        int width = vsapi->getFrameWidth(dst, 0);
        int height = vsapi->getFrameHeight(dst, 0);
        int hblockx = blockx / 2;

        for (int y = 0; y < height; y++) {
            int ydest = y / blocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx) {
                int m = std::min(width, x + hblockx);

                if (y % blocky == 0 || y % blocky == blocky - 1) {
                    for (int xl = x; xl < m; xl++) {
                        if (mask[ydest * nxblocks + xdest] == 1) {
                            dstp[xl] = markLow;
                        } else if (mask[ydest * nxblocks + xdest] == 2) {
                            dstp[xl] = markHigh;
                        }
                    }
                }

                if (mask[ydest * nxblocks + xdest] == 1) {
                    dstp[x] = markLow;
                    dstp[m - 1] = markLow;
                } else if (mask[ydest * nxblocks + xdest] == 2) {
                    dstp[x] = markHigh;
                    dstp[m - 1] = markHigh;
                }

                xdest++;
            }

            dstp += stride;
        }
    }
}

typedef struct {
    VSNode *node;
    const VSVideoInfo *vi;
    int order;
    int offset;
    int64_t dupthresh;
    int tratio;
    bool show;
} DotKillTData;

template<typename T>
static void dotKillTProcess(VSFrame *outframe, const VSFrame *srcpp, const VSFrame *srcp, const VSFrame *srcc, const VSFrame *srcn, const VSFrame *srcnn,
                            const DotKillTData *d, int pattern, int nxblocks, int nyblocks, int shift, const VSAPI *vsapi) {
    if (pattern == 0) {
        // 1
        applyDotcrawInverse<T>(srcc, srcn, outframe, d->order, shift, vsapi);

        // 2
        std::vector<int64_t> maskprev1(nxblocks * nyblocks);
        std::vector<int64_t> masknext1(nxblocks * nyblocks);
        std::vector<uint8_t> mask(nxblocks * nyblocks);
        calcDiffMetric<T>(srcp, srcn, maskprev1.data(), nxblocks, d->order, vsapi);
        calcDiffMetric<T>(srcc, srcnn, masknext1.data(), nxblocks, d->order, vsapi);

        diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

        applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, !d->order, d->show, shift, vsapi);
    } else if (pattern == 1) {
        // 1
        applyFieldBlend<T>(srcc, srcp, outframe, d->order, vsapi);

        // 2
        std::vector<int64_t> maskprev1(nxblocks * nyblocks);
        std::vector<int64_t> masknext1(nxblocks * nyblocks);
        std::vector<uint8_t> mask(nxblocks * nyblocks);
        calcDiffMetric<T>(srcp, srcn, maskprev1.data(), nxblocks, d->order, vsapi);
        calcDiffMetric<T>(srcc, srcnn, masknext1.data(), nxblocks, !d->order, vsapi);

        diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

        applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, !d->order, d->show, shift, vsapi);
    } else if (pattern == 2) {
        // 2
        // done before the masking below even though it is the second field, because
        // the two write disjoint fields from source frames only while the show
        // overlay spans every row and has to be drawn last to survive
        applyFieldBlend<T>(srcc, srcn, outframe, !d->order, vsapi);

        // 1
        std::vector<int64_t> maskprev1(nxblocks * nyblocks);
        std::vector<int64_t> masknext1(nxblocks * nyblocks);
        std::vector<uint8_t> mask(nxblocks * nyblocks);
        calcDiffMetric<T>(srcc, srcpp, maskprev1.data(), nxblocks, d->order, vsapi);
        calcDiffMetric<T>(srcp, srcn, masknext1.data(), nxblocks, !d->order, vsapi);

        diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

        applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, d->order, d->show, shift, vsapi);
    } else if (pattern == 3) {
        // 2
        applyDotcrawInverse<T>(srcc, srcp, outframe, !d->order, shift, vsapi);

        // 1
        std::vector<int64_t> maskprev1(nxblocks * nyblocks);
        std::vector<int64_t> masknext1(nxblocks * nyblocks);
        std::vector<uint8_t> mask(nxblocks * nyblocks);
        calcDiffMetric<T>(srcc, srcpp, maskprev1.data(), nxblocks, !d->order, vsapi);
        calcDiffMetric<T>(srcp, srcn, masknext1.data(), nxblocks, !d->order, vsapi);

        diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

        applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, d->order, d->show, shift, vsapi);
    } else if (pattern == 4) {
        // 1
        {
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric<T>(srcc, srcpp, maskprev1.data(), nxblocks, !d->order, vsapi);
            calcDiffMetric<T>(srcc, srcnn, masknext1.data(), nxblocks, d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

            applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, d->order, d->show, shift, vsapi);
        }

        // 2
        {
            std::vector<int64_t> maskprev1(nxblocks * nyblocks);
            std::vector<int64_t> masknext1(nxblocks * nyblocks);
            std::vector<uint8_t> mask(nxblocks * nyblocks);
            calcDiffMetric<T>(srcpp, srcc, maskprev1.data(), nxblocks, !d->order, vsapi);
            calcDiffMetric<T>(srcc, srcnn, masknext1.data(), nxblocks, d->order, vsapi);

            diffMetricToMask(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, d->dupthresh, d->tratio);

            applyTemporalMask<T>(outframe, srcc, srcp, srcn, mask.data(), nxblocks, !d->order, d->show, shift, vsapi);
        }
    }
}

static const VSFrame *VS_CC dotKillTGetFrame(int n, int activationReason, void *instanceData, [[maybe_unused]] void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    DotKillTData *d = reinterpret_cast<DotKillTData*>(instanceData);

    // a block covers blockx/2 columns and blocky frame rows, the latter being blocky/2 rows of a single field
    int nxblocks = (d->vi->width + blockx / 2 - 1) / (blockx / 2);
    int nyblocks = (d->vi->height + blocky - 1) / blocky;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(std::max(n - 2, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
        vsapi->requestFrameFilter(n + 2, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *srcpp = vsapi->getFrameFilter(std::max(n - 2, 0), d->node, frameCtx);
        const VSFrame *srcp = vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        const VSFrame *srcc = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *srcn = vsapi->getFrameFilter(n + 1, d->node, frameCtx);
        const VSFrame *srcnn = vsapi->getFrameFilter(n + 2, d->node, frameCtx);

        // first two fields are duplicates, meaning that offset 0, 1 and 5, 6 are used in the various calculations
        // fields 2, 3 and 4 needs to have determined how many blocks are consecutively static
        // note that comparisons are only run on a single since we can in most comparisons can eliminate
        // the dotcrawl from the equation

        // 0-1 2 3 4 5-6

        // the complementary field can likewise be used for movement detection

        /*
            FIELD OFFSETS
             +  -  +  -  +  -
            -1  0  1  2  3  4
            A1 B1 B1 C1 D1 E1
            A2 B2 C2 D2 D2 E2
        */

        VSFrame *outframe = vsapi->copyFrame(srcc, core);

        vsapi->mapSetInt(vsapi->getFramePropertiesRW(outframe), "DotKillTOffset", (n + d->offset) % 5, maReplace);

        int shift = d->vi->format.bitsPerSample - 8;

        if (d->vi->format.bytesPerSample == 1)
            dotKillTProcess<uint8_t>(outframe, srcpp, srcp, srcc, srcn, srcnn, d, (n + d->offset) % 5, nxblocks, nyblocks, shift, vsapi);
        else
            dotKillTProcess<uint16_t>(outframe, srcpp, srcp, srcc, srcn, srcnn, d, (n + d->offset) % 5, nxblocks, nyblocks, shift, vsapi);

        vsapi->freeFrame(srcpp);
        vsapi->freeFrame(srcp);
        vsapi->freeFrame(srcc);
        vsapi->freeFrame(srcn);
        vsapi->freeFrame(srcnn);

        return outframe;
    }

    return nullptr;
}

static void VS_CC dotKillTFree(void *instanceData, [[maybe_unused]] VSCore *core, const VSAPI *vsapi) {
    DotKillTData *d = reinterpret_cast<DotKillTData*>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC dotKillTCreate(const VSMap *in, VSMap *out, [[maybe_unused]] void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<DotKillTData> d(new DotKillTData());

    int err;
    d->node = vsapi->mapGetNode(in, "clip", 0, nullptr);
    d->vi = vsapi->getVideoInfo(d->node);
    d->offset = vsapi->mapGetIntSaturated(in, "offset", 0, &err);
    d->order = !!vsapi->mapGetInt(in, "order", 0, &err);
    int dupthresh = vsapi->mapGetIntSaturated(in, "dupthresh", 0, &err);
    if (err || dupthresh < 0)
        dupthresh = 64;
    // an 8 bit difference can never exceed 255, the headroom above that only exists
    // so the squaring below cannot overflow for absurd arguments
    dupthresh = std::min(dupthresh, 65535);
    d->tratio = vsapi->mapGetIntSaturated(in, "tratio", 0, &err);
    if (err || d->tratio < 1)
        d->tratio = 3;
    d->show = !!vsapi->mapGetInt(in, "show", 0, &err);

    if (!isSupportedFormat(d->vi, false)) {
        vsapi->mapSetError(out, "DotKillT: only constant dimension 8-16 bit integer YUV420 supported");
        vsapi->freeNode(d->node);
        return;
    }

    // dupthresh is given in 8 bit units so it has to be scaled to the sample range
    // before squaring, otherwise deeper input would never register as a duplicate
    int64_t scaledthresh = static_cast<int64_t>(dupthresh) << (d->vi->format.bitsPerSample - 8);
    d->dupthresh = scaledthresh * scaledthresh;

    if (d->offset < 0 || d->offset > 4) {
        vsapi->mapSetError(out, "DotKillT: offset must be between 0 and 4");
        vsapi->freeNode(d->node);
        return;
    }

    // diffMetricToMask needs at least one interior block plus a block of edge extension on every side
    if (d->vi->width < 3 * (blockx / 2) || d->vi->height < 3 * blocky) {
        vsapi->mapSetError(out, "DotKillT: clip dimensions must be at least 24x24");
        vsapi->freeNode(d->node);
        return;
    }

    VSFilterDependency deps[] = { {d->node, rpGeneral} };
    vsapi->createVideoFilter(out, "DotKillT", d->vi, dotKillTGetFrame, dotKillTFree, fmParallel, deps, 1, d.get(), core);
    d.release();
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("com.vapoursynth.dotkill", "dotkill", "VapourSynth DotKill", VS_MAKE_VERSION(4, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("DotKillS", "clip:vnode;iterations:int:opt;", "clip:vnode", dotKillSCreate, 0, plugin);
    vspapi->registerFunction("DotKillZ", "clip:vnode;order:int:opt;offset:int:opt;", "clip:vnode", dotKillZCreate, 0, plugin);
    vspapi->registerFunction("DotKillT", "clip:vnode;order:int:opt;offset:int:opt;dupthresh:int:opt;tratio:int:opt;show:int:opt;", "clip:vnode", dotKillTCreate, 0, plugin);
}
