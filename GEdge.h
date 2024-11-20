#include "include/GBitmap.h"
#include "include/GMath.h"
#include "include/GPoint.h"
struct GEdge {
    // X = mY + b
    float m;
    float b;
    GPoint top;
    GPoint bottom;
    int winding;

    bool horizontal() { return GRoundToInt(top.y) == GRoundToInt(bottom.y); }

    bool valid(float y) const { return y >= GRoundToInt(top.y) && y <= GRoundToInt(bottom.y); }

    void update() {
        if (top.y > bottom.y) {
            std::swap(top, bottom);
        }
        m = (top.x - bottom.x) / (top.y - bottom.y);
        b = top.x - m * top.y;
    }

    float getX(float y) { return m * y + b; }
    float getY(float x) { return (x - b) / m; }
};

inline GEdge createEdge(GPoint p1, GPoint p2) {
    GPoint top, bottom;
    int winding;
    if (GRoundToInt(p1.y) < GRoundToInt(p2.y)) {
        top = p1;
        bottom = p2;
        winding = -1;
    } else {
        top = p2;
        bottom = p1;
        winding = 1;
    }

    auto m = (top.x - bottom.x) / (top.y - bottom.y);
    auto b = top.x - m * top.y;
    auto edge = GEdge{m, b, top, bottom, winding};
    return edge;
}

inline void clipEdgeTo(std::vector<GEdge>& edges, const GBitmap bitmap, GEdge edge) {
    auto minX = 0;
    auto maxX = bitmap.width() - 1;
    auto minY = 0;
    auto maxY = bitmap.height() - 1;
    if (GRoundToInt(edge.top.y) < minY && GRoundToInt(edge.bottom.y) < minY) {
        return;
    }
    if (GRoundToInt(edge.top.y) > maxY && GRoundToInt(edge.bottom.y) > maxY) {
        return;
    }
    if (GRoundToInt(edge.top.y) < minY) {
        edge.top.x = edge.top.x + (edge.m * (minY - edge.top.y));
        edge.top.y = minY;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.bottom.y) > maxY) {
        edge.bottom.x = edge.bottom.x + (edge.m * (maxY - edge.bottom.y));
        edge.bottom.y = maxY;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.top.x) < minX && GRoundToInt(edge.bottom.x) < minX) {
        edge.top.x = minX;
        edge.bottom.x = minX;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.top.x) > maxX && GRoundToInt(edge.bottom.x) > maxX) {
        edge.top.x = maxX;
        edge.bottom.x = maxX;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.top.x) < minX) {
        auto newEdge = createEdge({static_cast<float>(minX), edge.top.y},
                                  {static_cast<float>(minX), (minX - edge.b) / edge.m});
        clipEdgeTo(edges, bitmap, newEdge);
        edge.top.x = minX;
        edge.top.y = (minX - edge.b) / edge.m;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.bottom.x) < minX) {
        auto newEdge = createEdge({static_cast<float>(minX), edge.bottom.y},
                                  {static_cast<float>(minX), (minX - edge.b) / edge.m});
        clipEdgeTo(edges, bitmap, newEdge);
        edge.bottom.x = minX;
        edge.bottom.y = (minX - edge.b) / edge.m;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.top.x) > maxX) {
        auto newEdge = createEdge({static_cast<float>(maxX), edge.top.y},
                                  {static_cast<float>(maxX), (maxX - edge.b) / edge.m});
        clipEdgeTo(edges, bitmap, newEdge);
        edge.top.x = maxX;
        edge.top.y = (maxX - edge.b) / edge.m;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    if (GRoundToInt(edge.bottom.x) > maxX) {
        auto newEdge = createEdge({static_cast<float>(maxX), edge.bottom.y},
                                  {static_cast<float>(maxX), (maxX - edge.b) / edge.m});
        clipEdgeTo(edges, bitmap, newEdge);
        edge.bottom.x = maxX;
        edge.bottom.y = (maxX - edge.b) / edge.m;
        clipEdgeTo(edges, bitmap, edge);
        return;
    }
    edge.update();
    if (!edge.horizontal()) {
        edges.push_back(edge);
    }
}
