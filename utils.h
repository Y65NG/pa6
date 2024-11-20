#include "GEdge.h"
#include "include/GBitmap.h"
#include "include/GBlendMode.h"
#include "include/GColor.h"
#include "include/GMath.h"
#include "include/GMatrix.h"
#include "include/GPaint.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GShader.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

/** Helper Functions */
inline unsigned int div255(unsigned int n) {
    return ((n + 128) * 257) >> 16;
    // return (n * ((1 << 16) + (1 << 8) + 1) + (1 << 23)) >> 24;
}

inline GPixel colorToPixel(GColor color) {
    color.a = std::max(0.0f, std::min(1.0f, color.a));
    color.r = std::max(0.0f, std::min(1.0f, color.r));
    color.g = std::max(0.0f, std::min(1.0f, color.g));
    color.b = std::max(0.0f, std::min(1.0f, color.b));
    return GPixel_PackARGB(GRoundToInt(color.a * 255), GRoundToInt(color.r * color.a * 255),
                           GRoundToInt(color.g * color.a * 255),
                           GRoundToInt(color.b * color.a * 255));
}

inline void src_row(GPixel* row_ptr, int count, GColor color) {
    auto pixelToFill = colorToPixel(color);
    for (auto x = 0; x < count; x++) {
        row_ptr[x] = pixelToFill;
    }
}

inline void blendClear(GPixel src, GPixel* dst) { *dst = 0; }

inline void blendSrc(GPixel src, GPixel* dst) { *dst = src; }

inline void blendDst(GPixel src, GPixel* dst) {}

inline void blendSrcOver(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(GPixel_GetA(src) + div255((255 - GPixel_GetA(src)) * GPixel_GetA(*dst)),
                           GPixel_GetR(src) + div255((255 - GPixel_GetA(src)) * GPixel_GetR(*dst)),
                           GPixel_GetG(src) + div255((255 - GPixel_GetA(src)) * GPixel_GetG(*dst)),
                           GPixel_GetB(src) + div255((255 - GPixel_GetA(src)) * GPixel_GetB(*dst)));
}

inline void blendDstOver(GPixel src, GPixel* dst) {
    *dst =
        GPixel_PackARGB(GPixel_GetA(*dst) + div255((255 - GPixel_GetA(*dst)) * GPixel_GetA(src)),
                        GPixel_GetR(*dst) + div255((255 - GPixel_GetA(*dst)) * GPixel_GetR(src)),
                        GPixel_GetG(*dst) + div255((255 - GPixel_GetA(*dst)) * GPixel_GetG(src)),
                        GPixel_GetB(*dst) + div255((255 - GPixel_GetA(*dst)) * GPixel_GetB(src)));
}

inline void blendSrcIn(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(
        div255(GPixel_GetA(*dst) * GPixel_GetA(src)), div255(GPixel_GetA(*dst) * GPixel_GetR(src)),
        div255(GPixel_GetA(*dst) * GPixel_GetG(src)), div255(GPixel_GetA(*dst) * GPixel_GetB(src)));
}

inline void blendDstIn(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(
        div255(GPixel_GetA(src) * GPixel_GetA(*dst)), div255(GPixel_GetA(src) * GPixel_GetR(*dst)),
        div255(GPixel_GetA(src) * GPixel_GetG(*dst)), div255(GPixel_GetA(src) * GPixel_GetB(*dst)));
}

inline void blendSrcOut(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(div255((255 - GPixel_GetA(*dst)) * GPixel_GetA(src)),
                           div255((255 - GPixel_GetA(*dst)) * GPixel_GetR(src)),
                           div255((255 - GPixel_GetA(*dst)) * GPixel_GetG(src)),
                           div255((255 - GPixel_GetA(*dst)) * GPixel_GetB(src)));
}

inline void blendDstOut(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(div255((255 - GPixel_GetA(src)) * GPixel_GetA(*dst)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetR(*dst)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetG(*dst)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetB(*dst)));
}

inline void blendSrcATop(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(
        div255(GPixel_GetA(*dst) * GPixel_GetA(src) + (255 - GPixel_GetA(src)) * GPixel_GetA(*dst)),
        div255(GPixel_GetA(*dst) * GPixel_GetR(src) + (255 - GPixel_GetA(src)) * GPixel_GetR(*dst)),
        div255(GPixel_GetA(*dst) * GPixel_GetG(src) + (255 - GPixel_GetA(src)) * GPixel_GetG(*dst)),
        div255(GPixel_GetA(*dst) * GPixel_GetB(src) +
               (255 - GPixel_GetA(src)) * GPixel_GetB(*dst)));
}

inline void blendDstATop(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(
        div255(GPixel_GetA(src) * GPixel_GetA(*dst) + (255 - GPixel_GetA(*dst)) * GPixel_GetA(src)),
        div255(GPixel_GetA(src) * GPixel_GetR(*dst) + (255 - GPixel_GetA(*dst)) * GPixel_GetR(src)),
        div255(GPixel_GetA(src) * GPixel_GetG(*dst) + (255 - GPixel_GetA(*dst)) * GPixel_GetG(src)),
        div255(GPixel_GetA(src) * GPixel_GetB(*dst) +
               (255 - GPixel_GetA(*dst)) * GPixel_GetB(src)));
}

inline void blendXor(GPixel src, GPixel* dst) {
    *dst = GPixel_PackARGB(div255((255 - GPixel_GetA(src)) * GPixel_GetA(*dst) +
                                  (255 - GPixel_GetA(*dst)) * GPixel_GetA(src)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetR(*dst) +
                                  (255 - GPixel_GetA(*dst)) * GPixel_GetR(src)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetG(*dst) +
                                  (255 - GPixel_GetA(*dst)) * GPixel_GetG(src)),
                           div255((255 - GPixel_GetA(src)) * GPixel_GetB(*dst) +
                                  (255 - GPixel_GetA(*dst)) * GPixel_GetB(src)));
}

// rewrite using function pointer
const auto blendFuncs = std::vector<std::function<void(GPixel, GPixel*)>>{
    blendClear, blendSrc,    blendDst,    blendSrcOver, blendDstOver, blendSrcIn,
    blendDstIn, blendSrcOut, blendDstOut, blendSrcATop, blendDstATop, blendXor};

inline void blend_row(int x, int y, int count, GPaint paint, GBitmap fDevice, GMatrix ctm) {
    assert(count >= 0);
    auto row_ptr = fDevice.getAddr(x, y);
    auto mode = paint.getBlendMode();
    if (paint.peekShader()) {
        // std::cout << "shader" << std::endl;
        auto shader = paint.shareShader();
        // shader->setContext(ctm);
        GPixel pixelsToFill[count];
        shader->shadeRow(x, y, count, pixelsToFill);
        if (shader->isOpaque()) {
            if (mode == GBlendMode::kSrcOver) {
                mode = GBlendMode::kSrc;
            }
            if (mode == GBlendMode::kDstIn) {
                mode = GBlendMode::kDst;
            }
            if (mode == GBlendMode::kDstOut) {
                mode = GBlendMode::kClear;
            }
            if (mode == GBlendMode::kSrcATop) {
                mode = GBlendMode::kSrcIn;
            }
            if (mode == GBlendMode::kDstATop) {
                mode = GBlendMode::kDstOver;
            }
        }
        auto modeIdx = static_cast<int>(mode);

        for (auto i = 0; i < count; i++) {
            blendFuncs[modeIdx](pixelsToFill[i], &row_ptr[i]);
        }
    } else {
        // std::cout << "no shader" << std::endl;
        auto color = paint.getColor();
        auto pixelToFill = colorToPixel(color);
        if (color.a == 1) {
            switch (mode) {
            case GBlendMode::kSrcOver:
                mode = GBlendMode::kSrc;
                break;
            case GBlendMode::kDstIn:
                mode = GBlendMode::kDst;
                break;
            case GBlendMode::kDstOut:
                mode = GBlendMode::kClear;
                break;
            case GBlendMode::kSrcATop:
                mode = GBlendMode::kSrcIn;
                break;
            case GBlendMode::kDstATop:
                mode = GBlendMode::kDstOver;
                break;
            default:
                break;
            }

        } else if (color.a == 0) {
            switch (mode) {
            case GBlendMode::kSrc:
            case GBlendMode::kSrcIn:
            case GBlendMode::kDstIn:
            case GBlendMode::kSrcOut:
            case GBlendMode::kDstATop:
                mode = GBlendMode::kClear;
                break;
            case GBlendMode::kSrcOver:
            case GBlendMode::kDstOver:
            case GBlendMode::kDstOut:
            case GBlendMode::kSrcATop:
            case GBlendMode::kXor:
                mode = GBlendMode::kDst;
                break;
            default:
                break;
            }
        }
        auto modeIdx = static_cast<int>(mode);
        for (auto i = 0; i < count; i++) {
            blendFuncs[modeIdx](pixelToFill, &row_ptr[i]);
        }
    }
}

inline GPoint lerpPoint(const GPoint verts[4], float u, float v) {
    // auto p0 = verts[0] * (1 - u) + verts[1] * u;
    // auto p1 = verts[3] * (1 - u) + verts[2] * u;
    // return p0 * (1 - v) + p1 * v;
    return (1 - u) * (1 - v) * verts[0] + u * (1 - v) * verts[1] + (1 - u) * v * verts[3] +
           u * v * verts[2];
}

inline GColor lerpColor(const GColor colors[4], float u, float v) {
    // auto c0 = colors[0] * (1 - u) + colors[1] * u;
    // auto c1 = colors[3] * (1 - u) + colors[2] * u;
    // return c0 * (1 - v) + c1 * v;
    return (1 - u) * (1 - v) * colors[0] + u * (1 - v) * colors[1] + (1 - u) * v * colors[3] +
           u * v * colors[2];
}

class MyTrigGradientShader : public GShader {
public:
    MyTrigGradientShader(GPoint p0, GPoint p1, GPoint p2, GColor c0, GColor c1, GColor c2)
        : p0(p0), p1(p1), p2(p2), c0(c0), c1(c1), c2(c2) {
        auto u = p1 - p0;
        auto v = p2 - p0;
        auto m = GMatrix(u, v, p0);
        localMatrix = m.invert().value();
    }

    bool isOpaque() override {
        return GRoundToInt(c0.a) == 1 && GRoundToInt(c1.a) == 1 && GRoundToInt(c2.a) == 1;
    }

    bool setContext(const GMatrix& ctm) override {
        auto result = ctm.invert();
        if (result.has_value()) {
            localMatrix = localMatrix * *result;
            return true;
        }
        return false;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {

        GPoint p[] = {static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f};
        localMatrix.mapPoints(p, 1);
        auto dc1 = c1 - c0;
        auto dc2 = c2 - c0;
        auto c = p[0].x * dc1 + p[0].y * dc2 + c0;
        auto dc = dc1 * localMatrix[0] + dc2 * localMatrix[1];
        for (auto i = 0; i < count; i++) {
            row[i] = colorToPixel(c);
            c += dc;
        }
    }

private:
    GPoint p0, p1, p2;
    GColor c0, c1, c2;
    GMatrix localMatrix;
};

inline std::shared_ptr<GShader> GCreateTrigGradient(GPoint p0, GPoint p1, GPoint p2, GColor c0,
                                                    GColor c1, GColor c2) {

    return std::make_unique<MyTrigGradientShader>(MyTrigGradientShader(p0, p1, p2, c0, c1, c2));
}

class MyProxyShader : public GShader {
public:
    MyProxyShader(GShader* realShader, GMatrix extraTransform)
        : realShader(realShader), extraTransform(extraTransform) {}

    bool isOpaque() override { return realShader->isOpaque(); }

    bool setContext(const GMatrix& ctm) override {
        return realShader->setContext(ctm * extraTransform);
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        realShader->shadeRow(x, y, count, row);
    }

private:
    GShader* realShader;
    GMatrix extraTransform;
};

inline std::shared_ptr<GShader> GCreateProxyShader(GShader* realShader, GMatrix extraTransform) {

    return std::make_unique<MyProxyShader>(MyProxyShader(realShader, extraTransform));
}

class MyCombineTrigShader : public GShader {
public:
    MyCombineTrigShader(GShader* trigGradientShader, GShader* proxyShader)
        : trigGradientShader(trigGradientShader), proxyShader(proxyShader) {}

    bool isOpaque() override { return trigGradientShader->isOpaque() && proxyShader->isOpaque(); }

    bool setContext(const GMatrix& ctm) override {

        return proxyShader->setContext(ctm) && trigGradientShader->setContext(ctm);
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        GPixel trigGradientRow[count];
        GPixel proxyRow[count];
        trigGradientShader->shadeRow(x, y, count, trigGradientRow);
        proxyShader->shadeRow(x, y, count, proxyRow);
        for (auto i = 0; i < count; i++) {
            row[i] =
                GPixel_PackARGB(div255(GPixel_GetA(trigGradientRow[i]) * GPixel_GetA(proxyRow[i])),
                                div255(GPixel_GetR(trigGradientRow[i]) * GPixel_GetR(proxyRow[i])),
                                div255(GPixel_GetG(trigGradientRow[i]) * GPixel_GetG(proxyRow[i])),
                                div255(GPixel_GetB(trigGradientRow[i]) * GPixel_GetB(proxyRow[i])));
        }
    }

private:
    GShader* trigGradientShader;
    GShader* proxyShader;
};

inline std::shared_ptr<GShader> GCreateCombineTrigShader(GShader* trigGradientShader,
                                                         GShader* proxyShader) {

    return std::make_unique<MyCombineTrigShader>(
        MyCombineTrigShader(trigGradientShader, proxyShader));
}
