#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GColor.h"
#include "include/GMath.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include <cmath>
#include <memory>
#include <vector>

class MyShader : public GShader {
public:
    MyShader(GBitmap bitmap, GMatrix localMatrix, GTileMode tileMode)
        : bitmap(bitmap), localMatrix(*localMatrix.invert()), tileMode(tileMode) {}
    bool isOpaque() override { return bitmap.isOpaque(); }

    bool setContext(const GMatrix& ctm) override {
        auto result = ctm.invert();
        if (result.has_value()) {
            localMatrix = localMatrix * *result;
            return true;
        }
        return false;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        auto px = localMatrix[0] * (x + 0.5) + localMatrix[2] * (y + 0.5) + localMatrix[4];
        auto py = localMatrix[1] * (x + 0.5) + localMatrix[3] * (y + 0.5) + localMatrix[5];
        for (int i = 0; i < count; ++i) {
            auto ix = px;
            auto iy = py;
            // clamp
            switch (tileMode) {
            case GTileMode::kClamp:
                if (ix < 0)
                    ix = 0;
                if (ix >= bitmap.width())
                    ix = bitmap.width() - 1;
                if (iy < 0)
                    iy = 0;
                if (iy >= bitmap.height())
                    iy = bitmap.height() - 1;
                break;
            case GTileMode::kRepeat:
                if (ix < 0)
                    ix = bitmap.width() - std::abs(GRoundToInt(ix)) % (bitmap.width()) - 1;
                else
                    ix = std::abs(GRoundToInt(ix)) % bitmap.width();
                if (iy < 0)
                    iy = bitmap.height() - std::abs(GRoundToInt(iy)) % (bitmap.height()) - 1;
                else
                    iy = std::abs(GRoundToInt(iy)) % bitmap.height();
                break;
            case GTileMode::kMirror:
                ix = std::abs(GRoundToInt(ix) % (2 * bitmap.width()));
                iy = std::abs(GRoundToInt(iy) % (2 * bitmap.height()));
                if (ix >= bitmap.width())
                    ix = 2 * bitmap.width() - ix - 1;
                if (iy >= bitmap.height())
                    iy = 2 * bitmap.height() - iy - 1;
                break;
            }

            // floor
            row[i] = *bitmap.getAddr(GFloorToInt(ix), GFloorToInt(iy));

            px += localMatrix[0];
            py += localMatrix[1];
        }
    }

private:
    GBitmap bitmap;
    GMatrix localMatrix;
    GTileMode tileMode;
};

std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix,
                                             GTileMode tileMode) {
    auto bitmapCopy = bitmap;
    auto localMatrixCopy = localMatrix;
    return std::make_shared<MyShader>(MyShader(bitmapCopy, localMatrixCopy, tileMode));
}

class MyLinearGradientShader : public GShader {
    GPoint p0, p1;
    GMatrix localMatrix;
    std::vector<GColor> colors;

public:
    MyLinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count,
                           GTileMode tileMode)
        : p0(p0), p1(p1), tileMode(tileMode) {
        auto dx = p1.x - p0.x;
        auto dy = p1.y - p0.y;
        localMatrix = GMatrix({static_cast<float>(count - 1), 0},
                              {0, static_cast<float>(count - 1)}, {0, 0}) *
                      *GMatrix({dx, dy}, {-dy, dx}, {p0.x, p0.y}).invert();
        this->colors = std::vector<GColor>(colors, colors + count);
    }

    bool isOpaque() override {
        for (auto& color : colors) {
            if (color.a < 1)
                return false;
        }
        return true;
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
        auto px = localMatrix[0] * (x + 0.5) + localMatrix[2] * (y + 0.5) + localMatrix[4];
        for (int i = 0; i < count; ++i) {
            auto ix = px;

            // clamp
            switch (tileMode) {
            case GTileMode::kClamp:
                if (ix < 0)
                    ix = 0;
                if (ix >= colors.size() - 1)
                    ix = colors.size() - 1;
                break;
            case GTileMode::kRepeat: {
                if (ix < 0)
                    ix = colors.size() + ix;
                ix = std::fmod(std::abs(ix), colors.size() - 1);
                break;
            }
            case GTileMode::kMirror: {
                ix = std::fmod(std::abs(ix), 2 * colors.size() - 2);
                // ix = std::abs(GRoundToInt(ix) % static_cast<int>(2 * colors.size() - 2));
                if (ix >= colors.size() - 1)
                    ix = 2 * colors.size() - 2 - ix;
                break;
            }
            }

            auto t = ix - GFloorToInt(ix);
            auto c = (1 - t) * colors[GFloorToInt(ix)] + t * colors[GFloorToInt(ix) + 1];

            // floor
            row[i] = GPixel_PackARGB(GRoundToInt(c.a * 255), GRoundToInt(c.r * c.a * 255),
                                     GRoundToInt(c.g * c.a * 255), GRoundToInt(c.b * c.a * 255));

            px += localMatrix[0];
        }
    }

private:
    GTileMode tileMode;
};

std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[],
                                               int count, GTileMode tileMode) {

    return std::make_unique<MyLinearGradientShader>(
        MyLinearGradientShader(p0, p1, colors, count, tileMode));
}
