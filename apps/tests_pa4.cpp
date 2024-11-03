/**
 *  Copyright 2018 Mike Reed
 */

#include "../include/GPathBuilder.h"
#include "../include/GPoint.h"
#include "../include/GCanvas.h"
#include "tests.h"

static void test_path(GTestStats* stats) {
    GPoint pts[GPath::kMaxNextPoints];
    GPathBuilder bu;
    auto path = bu.detach();

    EXPECT_EQ(stats, path->bounds(), GRect::WH(0,0));
    EXPECT_FALSE(stats, GPath::Iter(*path).next(pts).has_value());

    bu.moveTo(10, 20);
    bu.lineTo(30, 0);
    path = bu.detach();
    EXPECT_EQ(stats, path->bounds(), GRect::LTRB(10, 0, 30, 20));

    GPath::Iter iter(*path);
    EXPECT_EQ(stats, iter.next(pts).value(), GPathVerb::kMove);
    EXPECT_EQ(stats, pts[0], (GPoint{10, 20}));
    EXPECT_EQ(stats, iter.next(pts).value(), GPathVerb::kLine);
    EXPECT_EQ(stats, pts[0], (GPoint{10, 20}));
    EXPECT_EQ(stats, pts[1], (GPoint{30, 0}));
    EXPECT_FALSE(stats, iter.next(pts).has_value());

    auto p2 = path->transform(GMatrix::Scale(0.5, 2));
    EXPECT_EQ(stats, p2->bounds(), GRect::LTRB(5, 0, 15, 40));

    const GPoint p[] = { {10, 10}, {20, 20}, {30, 30} };
    int N = 0;
    bu.addPolygon(p, 3);  N++;
    bu.addPolygon(p, 3);  N++;
    bu.addPolygon(p, 3);  N++;
    path = bu.detach();

    GPath::Edger edger(*path);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j) {
            auto v = edger.next(pts);
            EXPECT_TRUE(stats, v.has_value() &&
                               *v == GPathVerb::kLine &&
                               pts[0] == p[j] &&
                               pts[1] == p[(j+1)%3]);
        }
    }
    EXPECT_FALSE(stats, edger.next(pts).has_value());
}

static bool expect_iter(const GPath& path, const GPoint expected_pts[], int count) {
    GPath::Iter iter(path);
    GPoint pts[GPath::kMaxNextPoints];
    if (count > 0) {
        if (iter.next(pts) != GPathVerb::kMove || pts[0] != expected_pts[0]) {
            return false;
        }
        for (int i = 1; i < count; ++i) {
            if (iter.next(pts) != GPathVerb::kLine || pts[1] != expected_pts[i]) {
                return false;
            }
        }
    }
    return !iter.next(pts).has_value();
}

static void test_path_rect(GTestStats* stats) {
    GPathBuilder bu;
    auto p = bu.detach();
    EXPECT_TRUE(stats, GRect::WH(0, 0) == p->bounds());

    // allow 4 or 5 points in the rect-path (auto closed or manually closed)

    GRect r = GRect::LTRB(10, 20, 30, 40);
    bu.addRect(r, GPathDirection::kCW);
    p = bu.detach();

    EXPECT_TRUE(stats, p->bounds() == r);
    const GPoint pts0[] = {{10, 20}, {30, 20}, {30, 40}, {10, 40}, {10, 20}};
    EXPECT_TRUE(stats, expect_iter(*p, pts0, 4) ||
                       expect_iter(*p, pts0, 5));

    bu.addRect(r, GPathDirection::kCCW);
    p = bu.detach();

    EXPECT_TRUE(stats, p->bounds() == r);
    const GPoint pts1[] = {{10, 20}, {10, 40}, {30, 40}, {30, 20}, {10, 20}};
    EXPECT_TRUE(stats, expect_iter(*p, pts1, 4) ||
                       expect_iter(*p, pts1, 5));
}

static void test_path_poly(GTestStats* stats) {
    const GPoint pts[] = {{-10, -10}, {10, 0}, {0, 10}};
    GPathBuilder bu;
    bu.addPolygon(pts, GARRAY_COUNT(pts));
    auto p = bu.detach();
    EXPECT_TRUE(stats, GRect::LTRB(-10, -10, 10, 10) == p->bounds());
    EXPECT_TRUE(stats, expect_iter(*p, pts, GARRAY_COUNT(pts)));
}

static void test_path_transform(GTestStats* stats) {
    GPathBuilder bu;
    const GPoint pts[] = {{10, 20}, {30, 40}, {50, 60}, {70, 80}};
    GRect r = GRect::LTRB(10, 20, 70, 80);

    bu.addPolygon(pts, GARRAY_COUNT(pts));
    auto p = bu.detach();
    EXPECT_TRUE(stats, r == p->bounds());

    auto p2 = p->transform(GMatrix::Translate(-30, 20));
    assert(p.get() != p2.get());

    r = r.offset(-30, 20);
    EXPECT_TRUE(stats, r == p2->bounds());

    auto p3 = p->transform(GMatrix());
    assert(p.get() == p3.get());
}

static bool all_zeros(const void* data, size_t size) {
    const char* ptr = (const char*)data;
    for (size_t i = 0; i < size; ++i) {
        if (ptr[i] != 0) {
            return false;
        }
    }
    return true;
}

static void test_path_nodraw(GTestStats* stats) {
    const int w = 4, h = 4;
    GPixel pixels[w*h];
    GBitmap bm(w, h, w*4, pixels, false);
    auto canvas = GCreateCanvas(bm);

    const std::shared_ptr<GPath> nodrawpaths[] = {
        GPathBuilder().detach(),
        [](){ GPathBuilder p; p.moveTo(1,1); return p.detach(); }(),
        [](){ GPathBuilder p; p.moveTo(1,1); p.lineTo(2,2); return p.detach(); }(),
        [](){ GPathBuilder p; p.moveTo(1,1); p.lineTo(2,2); p.lineTo(1,1); return p.detach(); }(),
    };
    const GPaint paint;

    for (const auto& p : nodrawpaths) {
        memset(pixels, 0, sizeof(pixels));
        canvas->drawPath(*p, paint);
        EXPECT_TRUE(stats, all_zeros(pixels, sizeof(pixels)));
    }
}
