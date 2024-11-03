/**
 *  Copyright 2018 Mike Reed
 */

#include "../include/GPath.h"
#include "tests.h"

static void test_edger_quads(GTestStats* stats) {
    const GPoint p[] = { {10, 10}, {20, 20}, {30, 30} };
    int N = 0;

    GPathBuilder bu;
    bu.moveTo(p[0]); bu.quadTo(p[1], p[2]);  N++;
    bu.moveTo(p[0]); bu.quadTo(p[1], p[2]);  N++;
    bu.moveTo(p[0]); bu.quadTo(p[1], p[2]);  N++;
    auto path = bu.detach();

    GPath::Edger edger(*path);
    GPoint pts[3];
    for (int i = 0; i < N; ++i) {
        assert(edger.next(pts) == GPathVerb::kQuad);
        for (int j = 0; j < 3; ++j) {
            assert(pts[j] == p[j]);
        }
    }
    assert(edger.next(pts) == GPathVerb::kLine);
    assert(pts[0] == p[2]);
    assert(pts[1] == p[0]);
    assert(!edger.next(pts));
}

static void test_path_circle(GTestStats* stats) {
    GPathBuilder bu;
    
    for (GPathDirection dir : { GPathDirection::kCW, GPathDirection::kCCW }) {
        bu.reset();
        bu.addCircle({10, 10}, 10, dir);
        auto p = bu.detach();
        EXPECT_TRUE(stats, GRect::WH(20, 20) == p->bounds());
    }
}

static void test_path_transform2(GTestStats* stats) {
    GPathBuilder bu;

    bu.moveTo(10, 10); bu.lineTo(20, 10); bu.lineTo(20, 40); bu.quadTo(10, 40, 10, 10);
    auto p = bu.detach();
    EXPECT_TRUE(stats, p->bounds() == GRect::LTRB(10, 10, 20, 40));

    GMatrix mx = GMatrix::Scale(2, 3);

    auto q = p->transform(mx);

    GPath::Iter iter_p(*p),
                iter_q(*q);
    GPoint      pts_p[GPath::kMaxNextPoints],
                pts_q[GPath::kMaxNextPoints];

    for (;;) {
        auto np = iter_p.next(pts_p);
        auto nq = iter_q.next(pts_q);
        EXPECT_EQ(stats, np.has_value(), nq.has_value());
        if (!np.has_value() || !nq.has_value()) {
            break;
        }
        EXPECT_TRUE(stats, np.value() == nq.value());

        int count;
        switch (np.value()) {
            case GPathVerb::kMove:  count = 1; break;
            case GPathVerb::kLine:  count = 2; break;
            case GPathVerb::kQuad:  count = 3; break;
            case GPathVerb::kCubic: count = 4; break;
        }
        mx.mapPoints(pts_p, count);
        for (int i = 0; i < count; ++i) {
            EXPECT_TRUE(stats, pts_p[i] == pts_q[i]);
        }
    }
}

static bool nearly_eq(float a, float b) {
    return fabs(a - b) <= 0.0001f;
}

static bool nearly_eq(GPoint a, GPoint b) {
    return nearly_eq(a.x, b.x) && nearly_eq(a.y, b.y);
}

static bool nearly_eq(const GRect& a, const GRect& b) {
    return nearly_eq(a.left,   b.left)  &&
           nearly_eq(a.top,    b.top)   &&
           nearly_eq(a.right,  b.right) &&
           nearly_eq(a.bottom, b.bottom);
}

static void test_path_chop_quad(GTestStats* stats) {
    const GPoint src[] = {
        { 0, 100 }, { 0, 0 }, { 200, 0 },
    };
    GPoint dst[5];

    GPath::ChopQuadAt(src, dst, 0.5f);

    const GPoint expected[] = {
        { 0, 100 }, { 0, 50 }, { 50, 25 }, { 100, 0 }, { 200, 0 }
    };
    for (int i = 0; i < 5; ++i) {
        EXPECT_TRUE(stats, nearly_eq(dst[i], expected[i]));
    }
}

static void test_path_chop_cubic(GTestStats* stats) {
    const GPoint src[] = {
        { 0, 80 }, { 0, 0 }, { 160, 0 }, { 240, 80 },
    };
    GPoint dst[7];

    GPath::ChopCubicAt(src, dst, 0.25f);

    const GPoint expected[] = {
        { 0, 80 }, { 0, 60 }, { 10, 45 }, { 26.25, 35 },
        { 75, 5 }, { 180, 20 }, { 240, 80 },
    };
    for (int i = 0; i < 7; ++i) {
        EXPECT_TRUE(stats, nearly_eq(dst[i], expected[i]));
    }
}

static void test_path_bounds(GTestStats* stats) {
    auto check_bounds = [&](GRect bounds, GRect expected) {
        EXPECT_TRUE(stats, nearly_eq(bounds, expected));
    };

    GPathBuilder bu;
    check_bounds(bu.detach()->bounds(), {0,0,0,0});

    bu.moveTo(0, 0);
    bu.quadTo(150, 0, 75, 75);
    bu.quadTo(0, 150, 0, 0);
    check_bounds(bu.detach()->bounds(), {0, 0, 100, 100});

    bu.moveTo (50, 100);
    bu.cubicTo(100, 0, 0, 0, 50, 100);
    check_bounds(bu.detach()->bounds(), {35.5662f, 25, 64.4338f, 100});

    bu.moveTo (0, 50);
    bu.cubicTo(100, 0, 100, 100, 0, 50);
    check_bounds(bu.detach()->bounds(), {0, 35.5662f, 75, 64.4338f});
}
