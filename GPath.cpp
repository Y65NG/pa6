#include "include/GPath.h"
#include "GEdge.h"
#include "include/GMatrix.h"
#include "include/GPathBuilder.h"
#include "include/GPoint.h"
#include "include/GRect.h"
#include <algorithm>
#include <cmath>
#include <vector>

void update_bounds(GRect& bounds, const GPoint& p) {
    if (std::abs(p.x) == INFINITY || std::abs(p.y) == INFINITY) {
        return;
    }
    bounds.left = std::min(bounds.left, p.x);
    bounds.top = std::min(bounds.top, p.y);
    bounds.right = std::max(bounds.right, p.x);
    bounds.bottom = std::max(bounds.bottom, p.y);
}

GPoint lerp(const GPoint& p1, const GPoint& p2, float t) { return p1 + t * (p2 - p1); }

GPoint getQuadPoint(const GPoint src[3], float t) {
    auto A = src[0];
    auto B = src[1];
    auto C = src[2];
    auto D = lerp(A, B, t);
    auto E = lerp(B, C, t);
    auto P = lerp(D, E, t);
    return P;
}

GPoint getCubicPoint(const GPoint src[4], float t) {
    auto A = src[0];
    auto B = src[1];
    auto C = src[2];
    auto D = src[3];
    auto E = lerp(A, B, t);
    auto F = lerp(B, C, t);
    auto G = lerp(C, D, t);
    auto H = lerp(E, F, t);
    auto I = lerp(F, G, t);
    auto J = lerp(H, I, t);
    return J;
}

GRect GPath::bounds() const {
    if (fPts.empty()) {
        return GRect::WH(0, 0);
    }

    Edger edger(*this);
    GPoint pts[kMaxNextPoints];

    GRect bounds = GRect::LTRB(INFINITY, INFINITY, -INFINITY, -INFINITY);

    for (auto v = edger.next(pts); v.has_value(); v = edger.next(pts)) {
        switch (v.value()) {
        case kMove:
            update_bounds(bounds, pts[0]);
            break;
        case kLine:
            update_bounds(bounds, pts[0]);
            update_bounds(bounds, pts[1]);
            break;
        case kQuad: {
            auto A = pts[0];
            auto B = pts[1];
            auto C = pts[2];
            update_bounds(bounds, A);
            update_bounds(bounds, C);

            auto tx = (A.x - B.x) / (A.x - 2 * B.x + C.x); // t when dx = 0
            auto ty = (A.y - B.y) / (A.y - 2 * B.y + C.y); // t when dy = 0

            update_bounds(bounds, getQuadPoint(pts, tx));
            update_bounds(bounds, getQuadPoint(pts, ty));
            break;
        }
        case kCubic: {
            auto A = pts[0];
            auto B = pts[1];
            auto C = pts[2];
            auto D = pts[3];
            update_bounds(bounds, A);
            update_bounds(bounds, D);

            auto ax = -3 * A.x + 9 * B.x - 9 * C.x + 3 * D.x;
            auto bx = 6 * A.x - 12 * B.x + 6 * C.x;
            auto cx = -3 * A.x + 3 * B.x;
            auto qx = -0.5 * (bx + (bx >= 0 ? (std::sqrt(bx * bx - 4 * ax * cx))
                                            : (-std::sqrt(bx * bx - 4 * ax * cx))));
            auto tx1 = qx / ax;
            auto tx2 = cx / qx;

            update_bounds(bounds, getCubicPoint(pts, tx1));
            update_bounds(bounds, getCubicPoint(pts, tx2));

            auto ay = -3 * A.y + 9 * B.y - 9 * C.y + 3 * D.y;
            auto by = 6 * A.y - 12 * B.y + 6 * C.y;
            auto cy = -3 * A.y + 3 * B.y;
            auto qy = -0.5 * (by + (by >= 0 ? (std::sqrt(by * by - 4 * ay * cy))
                                            : (-std::sqrt(by * by - 4 * ay * cy))));
            auto ty1 = qy / ay;
            auto ty2 = cy / qy;

            update_bounds(bounds, getCubicPoint(pts, ty1));
            update_bounds(bounds, getCubicPoint(pts, ty2));

            break;
        }
        }
    }

    return bounds;
}

void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    auto A = src[0];
    auto B = src[1];
    auto C = src[2];
    auto D = lerp(A, B, t);
    auto E = lerp(B, C, t);
    auto P = lerp(D, E, t);
    dst[0] = A;
    dst[1] = D;
    dst[2] = P;
    dst[3] = E;
    dst[4] = C;
}

void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    auto A = src[0];
    auto B = src[1];
    auto C = src[2];
    auto D = src[3];
    auto E = lerp(A, B, t);
    auto F = lerp(B, C, t);
    auto G = lerp(C, D, t);
    auto H = lerp(E, F, t);
    auto I = lerp(F, G, t);
    auto J = lerp(H, I, t);
    dst[0] = A;
    dst[1] = E;
    dst[2] = H;
    dst[3] = J;
    dst[4] = I;
    dst[5] = G;
    dst[6] = D;
}

void GPathBuilder::addRect(const GRect& rect, GPathDirection direction) {
    moveTo(rect.x(), rect.y());
    if (direction == GPathDirection::kCW) {
        lineTo(rect.x() + rect.width(), rect.y());
        lineTo(rect.x() + rect.width(), rect.y() + rect.height());
        lineTo(rect.x(), rect.y() + rect.height());
    } else {
        lineTo(rect.x(), rect.y() + rect.height());
        lineTo(rect.x() + rect.width(), rect.y() + rect.height());
        lineTo(rect.x() + rect.width(), rect.y());
    }
}

void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    assert(count > 0);
    moveTo(pts[0]);
    for (int i = 1; i < count; ++i) {
        lineTo(pts[i]);
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection direction) {
    float k = 0.551915;
    auto mx = GMatrix::Translate(center.x, center.y) * GMatrix::Scale(radius, radius);
    if (direction == GPathDirection::kCW) {
        auto pts = std::vector<GPoint>{
            {1, 0},  {1, k},   {k, 1},   {0, 1},  {-k, 1}, {-1, k},
            {-1, 0}, {-1, -k}, {-k, -1}, {0, -1}, {k, -1}, {1, -k},
        };
        mx.mapPoints(pts.data(), pts.size());
        moveTo(pts[0]);
        cubicTo(pts[1], pts[2], pts[3]);
        cubicTo(pts[4], pts[5], pts[6]);
        cubicTo(pts[7], pts[8], pts[9]);
        cubicTo(pts[10], pts[11], pts[0]);

        // moveTo(1, 0);
        // cubicTo({1, k}, {k, 1}, {0, 1});
        // cubicTo({-k, 1}, {-1, k}, {-1, 0});
        // cubicTo({-1, -k}, {-k, -1}, {0, -1});
        // cubicTo({k, -1}, {1, -k}, {1, 0});
    } else {
        auto pts = std::vector<GPoint>{
            {0, 1},  {k, 1},   {1, k},   {1, 0},  {1, -k}, {k, -1},
            {0, -1}, {-k, -1}, {-1, -k}, {-1, 0}, {-1, k}, {-k, 1},
        };
        mx.mapPoints(pts.data(), pts.size());
        moveTo(pts[0]);
        cubicTo(pts[1], pts[2], pts[3]);
        cubicTo(pts[4], pts[5], pts[6]);
        cubicTo(pts[7], pts[8], pts[9]);
        cubicTo(pts[10], pts[11], pts[0]);
        // moveTo(0, 1);
        // cubicTo({k, 1}, {1, k}, {1, 0});
        // cubicTo({1, -k}, {k, -1}, {0, -1});
        // cubicTo({-k, -1}, {-1, -k}, {-1, 0});
        // cubicTo({-1, k}, {-k, 1}, {0, 1});
    }
    // mx.mapPoints(fPts.data(), fPts.data(), fPts.size());
    // transform(mx);
}
