/**
 *  Copyright 2015 Mike Reed
 */

#include "image.h"
#include "../include/GCanvas.h"
#include "../include/GBitmap.h"
#include "../include/GColor.h"
#include "../include/GMatrix.h"
#include "../include/GPathBuilder.h"
#include "../include/GPoint.h"
#include "../include/GRandom.h"
#include "../include/GRect.h"
#include <string>

static void make_star(GPathBuilder* path, int count, float anglePhase) {
    assert(count & 1);
    float da = 2 * gFloatPI * (count >> 1) / count;
    float angle = anglePhase;
    for (int i = 0; i < count; ++i) {
        GPoint p = { cosf(angle), sinf(angle) };
        i == 0 ? path->moveTo(p) : path->lineTo(p);
        angle += da;
    }
}

static void add_star(GPathBuilder* path, int count) {
    if (count & 1) {
        make_star(path, count, 0);
    } else {
        count >>= 1;
        make_star(path, count, 0);
        make_star(path, count, gFloatPI);
    }
}

static void make_diamonds(GPathBuilder* path) {
    const GPoint pts[] {
        { 0, 0.1f }, { -1, 5 }, { 0, 6 }, { 1, 5 },
    };
    float steps = 12;
    
    GMatrix matrix;
    for (float angle = 0; angle < 2*gFloatPI - 0.001f; angle += 2*gFloatPI/steps) {
        GPoint dst[4];
        GMatrix::Rotate(angle).mapPoints(dst, pts, 4);
        path->addPolygon(dst, 4);
    }
}

static void add_poly(GPathBuilder* bu, std::initializer_list<GPoint> list) {
    bu->addPolygon(list.begin(), list.size());
}

static void stars(GCanvas* canvas) {
    canvas->clear({1, 1, 1, 1});
    
    GMatrix scale = GMatrix::Scale(100, 100);
    
    GPathBuilder bu;
    add_star(&bu, 5);
    bu.transform(scale);
    canvas->translate(120, 120);
    canvas->drawPath(*bu.detach(), GPaint({1,0,0,1}));
    
    canvas->translate(250, 0);
    add_star(&bu, 6);
    bu.transform(scale);
    canvas->drawPath(*bu.detach(), GPaint({0,0,1,1}));
    
    canvas->translate(-250, 250);
    add_poly(&bu, {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}});
    add_poly(&bu, {{0, -1}, {-1,  }, {0, 1}, { 1, 0}});
    add_poly(&bu, {{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}});
    add_poly(&bu, {{0, -0.5}, {-0.5, 0}, {0, 0.5}, {0.5, 0}});

    bu.transform(scale);
    canvas->drawPath(*bu.detach(), GPaint());
    
    make_diamonds(&bu);
    bu.transform(GMatrix::Scale(20, 20));
    canvas->translate(250, 0);
    canvas->drawPath(*bu.detach(), GPaint({0, 1, 0, 1}));
}

struct ARGB {
    float a, r, g, b;

    operator GColor() const { return GColor{r, g, b, a}; }
};

static void draw_lion_bare(GCanvas* canvas) {
#include "lion.inc"
}

static void draw_lion(GCanvas* canvas) {
    canvas->translate(130, 40);
    canvas->scale(1.2f, 1.2f);
    draw_lion_bare(canvas);
}

static void draw_lion_head(GCanvas* canvas) {
    canvas->translate(-40, 10);
    canvas->scale(3, 3);
    draw_lion_bare(canvas);
}

static void draw_grad(GCanvas* canvas) {
    canvas->scale(2, 2);

    GRect r = GRect::XYWH(10, 15, 100, 75);
    const GColor colors[] = {
        {1,0,0,1}, {0,1,0,1}, {0,0,1,1},
    };
    auto sh = GCreateLinearGradient({r.left, r.top}, {r.right, r.bottom},
                                    colors, GARRAY_COUNT(colors));
    GPaint paint(sh);
    canvas->drawRect({-100, -100, 600, 600}, paint);
}

static void draw_mode_gradients(GCanvas* canvas, const GRect& bounds, GBlendMode mode) {
    outer_frame(canvas, bounds);

    auto dstsh = GCreateLinearGradient({bounds.left, bounds.top}, {bounds.left, bounds.bottom},
                                       {0, 0, 0, 0}, {1, 0, 0, 1});
    auto srcsh = GCreateLinearGradient({bounds.left, bounds.top}, {bounds.right, bounds.top},
                                       {0, 0, 0, 0}, {0, 0, 1, 1});
    GPaint paint;
    GPathBuilder bu;
    bu.addRect(bounds);

    paint.setBlendMode(GBlendMode::kSrc);
    paint.setShader(dstsh);
    canvas->drawRect(bounds, paint);

    paint.setBlendMode(mode);
    paint.setShader(srcsh);
    canvas->drawPath(*bu.detach(), paint);
}

static void draw_gradient_blendmodes(GCanvas* canvas) {
    draw_all_blendmodes(canvas, draw_mode_gradients);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

static void graph_path(GPathBuilder* path, int steps, float min, float max, float (*func)(float x)) {
    const float dx = (max - min) / (steps - 1);
    float x = min;
    path->moveTo(x, -func(x));
    for (int i = 1; i < steps; ++i) {
        x += dx;
        path->lineTo(x, -func(x));
    }
}

static void draw_graphs2(GCanvas* canvas) {
    GPathBuilder bu;

    float min = -gFloatPI;
    float max =  gFloatPI;
    canvas->save();
    canvas->translate(128, 128);
    canvas->scale(40, 60);
    graph_path(&bu, 30, min, max, [](float x) { return sinf(x); });
    auto sh = GCreateLinearGradient({min, 0}, {max, 0}, {0, 0, 1, 1}, {1, 0, 0, 1});
    canvas->drawPath(*bu.detach(), GPaint(sh));
    canvas->restore();

    GColor color = {0, 0.5f, 0, 1};
    min = -5*gFloatPI;
    max =  5*gFloatPI;
    canvas->save();
    canvas->translate(128, 40);
    canvas->scale(10, 40);
    graph_path(&bu, 70, min, max, [](float x) { return sinf(x)/x; });
    sh = GCreateLinearGradient({0, -1.f}, {0, 1.f}, &color, 1);
    canvas->drawPath(*bu.detach(), GPaint(sh));
    canvas->restore();

    min = 0;
    max = 1;
    canvas->save();
    canvas->translate(128, 250);
    canvas->scale(100, 100);
    graph_path(&bu, 40, min, max, [](float x) { return sqrtf(x); });
    bu.lineTo(1, 0);
    sh = GCreateLinearGradient({min,0}, {max,0}, {0.75f,0.75f,0.75f, 1}, {0,0,0,1});
    canvas->drawPath(*bu.detach(), GPaint(sh));
    canvas->restore();
}

