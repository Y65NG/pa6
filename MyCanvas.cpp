/*
 *  Copyright 2024 <me>
 */

#include "MyCanvas.h"
#include "include/GBitmap.h"
#include "include/GColor.h"
#include "include/GMath.h"
#include "include/GMatrix.h"
#include "include/GPaint.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "include/GPoint.h"
#include "include/GRect.h"
#include "include/GShader.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

void MyCanvas::clear(const GColor& color) {
    // your code here
    //* the order traversing the memory matters because of caching
    for (auto h = 0; h < fDevice.height(); h++) {
        auto pixel = fDevice.getAddr(0, h);
        src_row(pixel, fDevice.width(), color); //* speed up
    }
}

void MyCanvas::drawRect(const GRect& rect, const GPaint& paint) {
    auto roundedRect = rect.round();

    if (paint.peekShader() || ctm != GMatrix()) {
        auto points = std::vector<GPoint>{
            {static_cast<float>(roundedRect.left), static_cast<float>(roundedRect.top)},
            {static_cast<float>(roundedRect.right), static_cast<float>(roundedRect.top)},
            {static_cast<float>(roundedRect.right), static_cast<float>(roundedRect.bottom)},
            {static_cast<float>(roundedRect.left), static_cast<float>(roundedRect.bottom)}};
        drawConvexPolygon(points.data(), 4, paint);
        return;
    }
    auto l = std::max(0, roundedRect.left);
    auto r = std::min(fDevice.width(), roundedRect.right);
    if (l >= r)
        return;

    for (auto h = std::max(0, roundedRect.top); h < std::min(fDevice.height(), roundedRect.bottom);
         h++) {
        blend_row(l, h, r - l, paint, fDevice, ctm);
    }
}

void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {

    if (paint.peekShader()) {

        paint.peekShader()->setContext(ctm);
    }
    GPoint newPoints[count];
    auto edges = std::vector<GEdge>();
    if (ctm != GMatrix()) {
        ctm.mapPoints(newPoints, points, count);

        for (auto i = 0; i < count; i += 1) {
            auto edge = createEdge(newPoints[i], newPoints[i + 1 == count ? 0 : i + 1]);
            if (!edge.horizontal()) {
                clipEdgeTo(edges, fDevice, edge);
            }
        }
    } else {
        for (auto i = 0; i < count; i += 1) {
            auto edge = createEdge(points[i], points[i + 1 == count ? 0 : i + 1]);
            if (!edge.horizontal()) {
                clipEdgeTo(edges, fDevice, edge);
            }
        }
    }

    if (edges.size() < 2)
        return;

    std::sort(edges.begin(), edges.end(), [](GEdge a, GEdge b) { return a.top.y < b.top.y; });

    auto i = 0;
    auto j = 1;
    auto nextIdx = 2;

    for (int y = GRoundToInt(edges[0].top.y); y < GRoundToInt(edges[edges.size() - 1].bottom.y);
         ++y) {
        if (i < edges.size() && y >= edges[i].bottom.y) {
            i = nextIdx++;
        }
        if (j < edges.size() && y >= edges[j].bottom.y) {
            j = nextIdx++;
        }
        if (y < edges[i].top.y || y < edges[j].top.y)
            continue;
        if (i >= edges.size() || j >= edges.size())
            break;
        auto left = GRoundToInt(std::min(edges[i].getX(y + 0.5), edges[j].getX(y + 0.5)));
        // auto left = GFloorToInt(std::min(edges[i].getX(r), edges[j].getX(r)));
        auto right = GRoundToInt(std::max(edges[i].getX(y + 0.5), edges[j].getX(y + 0.5)));

        blend_row(left, y, right - left, paint, fDevice, ctm);
    }
    if (paint.peekShader()) {
        auto inv = ctm.invert();
        if (inv.has_value()) {
            paint.peekShader()->setContext(*ctm.invert());
        }
    }
}

void MyCanvas::save() { copies.push_back(GMatrix(ctm)); }

void MyCanvas::restore() {
    assert(!copies.empty());
    ctm = GMatrix(copies.back());
    copies.pop_back();
}

void MyCanvas::concat(const GMatrix& matrix) { ctm = ctm * matrix; }

void createQuadEdgesTo(std::vector<GEdge>& edges, const GPoint src[3], int numToChop,
                       const GBitmap fDevice) {
    if (numToChop == 0) {
        auto edge = createEdge(src[0], src[2]);
        if (!edge.horizontal()) {
            clipEdgeTo(edges, fDevice, edge);
        }
        return;
    }
    GPoint dst[5];
    GPath::ChopQuadAt(src, dst, 0.5);
    createQuadEdgesTo(edges, dst, numToChop - 1, fDevice);
    createQuadEdgesTo(edges, dst + 2, numToChop - 1, fDevice);
}

void createCubicEdgesTo(std::vector<GEdge>& edges, const GPoint src[4], int numToChop,
                        const GBitmap fDevice) {
    if (numToChop == 0) {
        auto edge = createEdge(src[0], src[3]);
        if (!edge.horizontal()) {
            clipEdgeTo(edges, fDevice, edge);
        }
        return;
    }
    GPoint dst[7];
    GPath::ChopCubicAt(src, dst, 0.5);
    createCubicEdgesTo(edges, dst, numToChop - 1, fDevice);
    createCubicEdgesTo(edges, dst + 3, numToChop - 1, fDevice);
}

void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {

    if (paint.peekShader()) {
        paint.peekShader()->setContext(ctm);
    }
    auto transformedPath = path.transform(ctm);
    GPath::Edger edger(*transformedPath);
    GPoint pts[GPath::kMaxNextPoints];
    auto edges = std::vector<GEdge>();
    while (auto v = edger.next(pts)) {
        switch (v.value()) {
        case GPathVerb::kLine: {
            auto edge = createEdge(pts[0], pts[1]);
            if (!edge.horizontal()) {
                clipEdgeTo(edges, fDevice, edge);
            }
            break;
        }
        case GPathVerb::kQuad: {
            auto A = pts[0];
            auto B = pts[1];
            auto C = pts[2];
            auto E = A - 2 * B + C;
            auto err = std::abs(E.length() / 4);

            int numSegs = GCeilToInt(std::sqrt(err * 4));
            int numToChop = GCeilToInt(std::log2(numSegs));
            createQuadEdgesTo(edges, pts, numToChop, fDevice);
            break;
        }
        case GPathVerb::kCubic: {
            auto A = pts[0];
            auto B = pts[1];
            auto C = pts[2];
            auto D = pts[3];
            auto E0 = A - 2 * B + C;
            auto E1 = B - 2 * C + D;
            GPoint E = {std::max(E0.x, E1.x), std::max(E0.y, E1.y)};
            auto err = std::abs(E.length());
            int numSegs = GCeilToInt(std::sqrt(3 * err));
            int numToChop = GCeilToInt(std::log2(numSegs));
            createCubicEdgesTo(edges, pts, numToChop, fDevice);
        }
        default:
            break;
        }
    }
    if (edges.size() < 2)
        return;

    std::sort(edges.begin(), edges.end(), [](GEdge a, GEdge b) {
        if (GRoundToInt(a.top.y) == GRoundToInt(b.top.y)) {
            return a.getX(a.top.y + 0.5) < b.getX(b.top.y + 0.5);
        }
        return GRoundToInt(a.top.y) < GRoundToInt(b.top.y);
    });

    int yUpper = GRoundToInt(edges[0].top.y);
    int yLower = 0;
    for (auto edge : edges) {
        yLower = std::max(yLower, GRoundToInt(edge.bottom.y));
    }
    for (int y = yUpper; y <= yLower + 1; ++y) {
        size_t i = 0;
        int w = 0;
        int L;
        float r = y + 0.5;
        while (i < edges.size() && edges[i].valid(r)) {
            int x = GRoundToInt(edges[i].getX(r));
            if (x < 0) {
                x = 0;
            }
            if (x >= fDevice.width()) {
                x = fDevice.width() - 1;
            }
            if (w == 0) {
                L = x;
            }
            w += edges[i].winding;
            if (w == 0) {
                blend_row(L, y, x - L, paint, fDevice, ctm);
            }

            if (edges[i].valid(r + 1)) {
                ++i;
            } else {
                edges.erase(edges.begin() + i);
            }
        }
        while (i < edges.size() && edges[i].valid(r + 1)) {
            ++i;
        }
        std::sort(edges.begin(), edges.begin() + i,
                  [r](auto a, auto b) { return a.getX(r + 1) < b.getX(r + 1); });
    }
    if (paint.peekShader()) {
        auto inv = ctm.invert();
        if (inv.has_value()) {
            paint.peekShader()->setContext(*ctm.invert());
        }
    }
}

void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count,
                        const int indices[], const GPaint& paint) {
    int n = 0;
    GPaint p = paint;
    for (int i = 0; i < count; ++i) {
        auto P0 = verts[indices[n + 0]];
        auto P1 = verts[indices[n + 1]];
        auto P2 = verts[indices[n + 2]];
        auto points = std::vector<GPoint>({P0, P1, P2});

        if (colors != nullptr && texs == nullptr) {
            // TrigGradient
            auto C0 = colors[indices[n + 0]];
            auto C1 = colors[indices[n + 1]];
            auto C2 = colors[indices[n + 2]];

            p.setShader(GCreateTrigGradient(P0, P1, P2, C0, C1, C2));

        } else if (colors == nullptr && texs != nullptr && paint.peekShader() != nullptr) {
            // MyProxyShader
            auto T0 = texs[indices[n + 0]];
            auto T1 = texs[indices[n + 1]];
            auto T2 = texs[indices[n + 2]];

            auto P = GMatrix(P1 - P0, P2 - P0, P0);
            auto T = GMatrix(T1 - T0, T2 - T0, T0);

            p.setShader(GCreateProxyShader(paint.peekShader(), P * T.invert().value()));
        } else if (colors != nullptr && texs != nullptr && paint.peekShader() != nullptr) {
            // CombineTrigShader
            auto C0 = colors[indices[n + 0]];
            auto C1 = colors[indices[n + 1]];
            auto C2 = colors[indices[n + 2]];

            auto T0 = texs[indices[n + 0]];
            auto T1 = texs[indices[n + 1]];
            auto T2 = texs[indices[n + 2]];

            auto P = GMatrix(P1 - P0, P2 - P0, P0);
            auto T = GMatrix(T1 - T0, T2 - T0, T0);

            // auto trigGradientShader = MyTrigGradientShader(P0, P1, P2, C0, C1, C2);
            // auto proxyShader = MyProxyShader(paint.peekShader(), P * T.invert().value());
            p.setShader(GCreateCombineTrigShader(GCreateTrigGradient(P0, P1, P2, C0, C1, C2), GCreateProxyShader(paint.peekShader(), P * T.invert().value())));
        }

        drawConvexPolygon(points.data(), 3, p);
        n += 3;
    }
}

// void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int
// count,
//                         const int indices[], const GPaint& paint) {
//     int n = 0;
//     GShader* originalShader = paint.peekShader();

//     for (int i = 0; i < count; ++i) {
//         GPoint trigPts[3] = {verts[indices[n + 0]], verts[indices[n + 1]], verts[indices[n +
//         2]]};

//         if (colors != nullptr && texs == nullptr) {
//             GColor triangleColors[3] = {colors[indices[n + 0]], colors[indices[n + 1]],
//                                         colors[indices[n + 2]]};
//             // pack
//             std::shared_ptr<GShader> shader =
//                 GCreateTrigGradient(trigPts[0], trigPts[1], trigPts[2], triangleColors[0],
//                                     triangleColors[1], triangleColors[2]);
//             drawConvexPolygon(trigPts, 3, GPaint(shader));
//         }

//         if (colors == nullptr && texs != nullptr && originalShader != nullptr) {
//             GPoint trigTexs[3] = {texs[indices[n + 0]], texs[indices[n + 1]], texs[indices[n +
//             2]]}; GMatrix P = {trigPts[1].x - trigPts[0].x, trigPts[2].x - trigPts[0].x,
//             trigPts[0].x,
//                          trigPts[1].y - trigPts[0].y, trigPts[2].y - trigPts[0].y, trigPts[0].y};
//             GMatrix T = {
//                 trigTexs[1].x - trigTexs[0].x, trigTexs[2].x - trigTexs[0].x, trigTexs[0].x,
//                 trigTexs[1].y - trigTexs[0].y, trigTexs[2].y - trigTexs[0].y, trigTexs[0].y};
//             GMatrix invTransform = P * T.invert().value();
//             // pack
//             std::shared_ptr<GShader> proxy =
//                 std::make_shared<MyProxyShader>(originalShader, invTransform);
//             drawConvexPolygon(trigPts, 3, GPaint(proxy));
//         }

//         if (colors != nullptr && texs != nullptr && originalShader != nullptr) {
//             GPoint trigTexs[3] = {texs[indices[n + 0]], texs[indices[n + 1]], texs[indices[n +
//             2]]}; GColor trigColors[3] = {colors[indices[n + 0]], colors[indices[n + 1]],
//                                     colors[indices[n + 2]]};
//             GMatrix P = {trigPts[1].x - trigPts[0].x, trigPts[2].x - trigPts[0].x, trigPts[0].x,
//                          trigPts[1].y - trigPts[0].y, trigPts[2].y - trigPts[0].y, trigPts[0].y};

//             GMatrix T = {
//                 trigTexs[1].x - trigTexs[0].x, trigTexs[2].x - trigTexs[0].x, trigTexs[0].x,
//                 trigTexs[1].y - trigTexs[0].y, trigTexs[2].y - trigTexs[0].y, trigTexs[0].y};
//             GMatrix invTransform = P * T.invert().value();
//             MyTrigGradientShader MyTrigGradientShader(trigPts[0], trigPts[1], trigPts[2],
//                                                       trigColors[0], trigColors[1],
//                                                       trigColors[2]);
//             MyProxyShader MyProxyShader(originalShader, invTransform);
//             // pack
//             std::shared_ptr<GShader> combined =
//                 std::make_shared<MyCombineTrigShader>(&MyTrigGradientShader, &MyProxyShader);
//             drawConvexPolygon(trigPts, 3, GPaint(combined));
//         }
//         n += 3;
//     }
// };

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                        int level, const GPaint& paint) {
    GPoint subVerts[4];
    GColor subColors[4];
    GPoint subTexs[4];

    int indices[6] = {0, 1, 3, 1, 2, 3};
    for (int u = 0; u < level + 1; ++u) {
        auto u0 = static_cast<float>(u) / (level + 1);
        auto u1 = static_cast<float>(u + 1) / (level + 1);
        for (int v = 0; v < level + 1; ++v) {
            auto v0 = static_cast<float>(v) / (level + 1);
            auto v1 = static_cast<float>(v + 1) / (level + 1);

            float uv[4][2] = {{u0, v0}, {u1, v0}, {u1, v1}, {u0, v1}};

            for (int i = 0; i < 4; ++i) {
                subVerts[i] = lerpPoint(verts, uv[i][0], uv[i][1]);
                if (colors != nullptr) {
                    subColors[i] = lerpColor(colors, uv[i][0], uv[i][1]);
                }
                if (texs != nullptr) {
                    subTexs[i] = lerpPoint(texs, uv[i][0], uv[i][1]);
                }
            }

            drawMesh(subVerts, colors != nullptr ? subColors : nullptr,
                     texs != nullptr ? subTexs : nullptr, 2, indices, paint);
        }
    }
}

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    canvas->translate(-40, 10);
    canvas->scale(3, 3);
    GPathBuilder bu;

    auto draw = [&](const GPoint pts[], size_t count, const GPaint& paint) {
        bu.addPolygon(pts, count);
        canvas->drawPath(*bu.detach(), paint);
    };
    GPaint paint({0, 0, 0, 1});
    const GPoint pts[] = {
        {41, 94}, {32, 110}, {23, 132}, {12, 163}, {6, 190},  {7, 217},  {5, 236}, {3, 247},
        {9, 230}, {12, 211}, {12, 185}, {18, 160}, {26, 134}, {35, 110}, {43, 99}, {41, 94},
    };
    draw(pts, GARRAY_COUNT(pts), paint);

    return "tears in rain";
}
