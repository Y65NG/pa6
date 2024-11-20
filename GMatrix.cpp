#include "include/GMatrix.h"
#include "include/GPoint.h"
#include "include/nonstd/optional.hpp"
#include <cmath>
#include <iostream>
GMatrix::GMatrix() : GMatrix(1, 0, 0, 0, 1, 0) {}

GMatrix GMatrix::Translate(float tx, float ty) { return GMatrix({1, 0}, {0, 1}, {tx, ty}); }

GMatrix GMatrix::Scale(float sx, float sy) { return GMatrix({sx, 0}, {0, sy}, {0, 0}); }

GMatrix GMatrix::Rotate(float radians) {
    return GMatrix({std::cos(radians), std::sin(radians)}, {-std::sin(radians), std::cos(radians)},
                   {0, 0});
}

GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {
    return GMatrix(a[0] * b[0] + a[2] * b[1], a[0] * b[2] + a[2] * b[3],
                   a[0] * b[4] + a[2] * b[5] + a[4], a[1] * b[0] + a[3] * b[1],
                   a[1] * b[2] + a[3] * b[3], a[1] * b[4] + a[3] * b[5] + a[5]);
}

// Calculate the inverse of 3 by 3 matrices where the last row is [0, 0, 1]
nonstd::optional<GMatrix> GMatrix::invert() const {
    auto a = fMat[0];
    auto b = fMat[1];
    auto c = fMat[2];
    auto d = fMat[3];
    auto e = fMat[4];
    auto f = fMat[5];
    auto det = a * d - c * b;
    if (det == 0)
        return {};

    return GMatrix(d / det, -c / det, (c * f - d * e) / det, -b / det, a / det,
                   (b * e - a * f) / det);
    // return GMatrix(d / det, -c / det, -e / det, -b / det, a / det, -f / det);
}

void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
    for (int i = 0; i < count; ++i) {
        GPoint newDst = GPoint{fMat[0] * src[i].x + fMat[2] * src[i].y + fMat[4],
                               fMat[1] * src[i].x + fMat[3] * src[i].y + fMat[5]};
        dst[i] = newDst;
    }
}
