/**
 *  Copyright 2018 Mike Reed
 */

#include "draw_curve.h"

template <typename DRAW> void spin(GCanvas* canvas, int N, DRAW draw) {
    for (int i = 0; i < N; ++i) {
        canvas->save();
        canvas->rotate(2 * M_PI * i / N);
        draw(canvas);
        canvas->restore();
    }
}

class CubicShape : public CurveShape {
    GPoint fPts[4];
public:
    CubicShape() {
        fPts[0] = { 20, 100 };
        fPts[1] = { 100, 20 };
        fPts[2] = { 160, 120 };
        fPts[3] = { 250, 50 };
        fColor = { 1, 0, 0, 1 };
    }

    void onDraw(GCanvas* canvas, const GPaint& paint) override {
        if (fShowPoints) {
            this->drawPoints(canvas, fPts, 4);
            return;
        }
        auto path = GPathBuilder::Build([this](GPathBuilder& bu) {
            bu.moveTo(fPts[0]);
            bu.cubicTo(fPts[1], fPts[2], fPts[3]);
        });

        if (fShowBounds) {
            canvas->drawRect(path->bounds(), GPaint({0,0,1,0.5f}));
        }
        canvas->drawPath(*path, paint);

        if (fShowHalves) {
            GPoint dst[7];
            GPath::ChopCubicAt(fPts, dst, fSubT);

            path = GPathBuilder::Build([dst](GPathBuilder& bu) {
                for (int i = 1; i <= 5; ++i) {
                    bu.addCircle(dst[i], 4);
                }
            });
            canvas->drawPath(*path, GPaint());
        }
    }

    void drawHilite(GCanvas* canvas) override {
        auto add_dashed_line = [](GPathBuilder* bu, GPoint a, GPoint b) {
            GVector v = b - a;
            const float len = v.length();
            const float dist = 16;
            const GVector dv = v * (dist / len);
            a = a + dv * 0.5f;
            const float R = 1.4f;
            for (float d = dist * 0.5f; d < len; d += dist) {
                bu->addRect(GRect::LTRB(a.x - R, a.y - R, a.x + R, a.y + R));
                a = a + dv;
            }
        };
        GPaint paint;
        GPathBuilder bu;
        const int count = GARRAY_COUNT(fPts);
        for (int i = 0; i < count; ++i) {
            bu.addCircle(fPts[i], 6);
            if (!fShowPoints && i < count - 1) {
                add_dashed_line(&bu, fPts[i], fPts[i+1]);
            }
        }
        if (fShowHalves) {
            GPoint dst[7];
            GPath::ChopCubicAt(fPts, dst, fSubT);
            add_dashed_line(&bu, dst[1], dst[2]);
            add_dashed_line(&bu, dst[2], dst[4]);
            add_dashed_line(&bu, dst[4], dst[5]);
        }
        canvas->drawPath(*bu.detach(), paint);

        this->Shape::drawShaderHilite(canvas);
    }

    GRect getRect() override {
        GPathBuilder bu;
        bu.moveTo(fPts[0]);
        bu.cubicTo(fPts[1], fPts[2], fPts[3]);
        return bu.detach()->bounds();
    }

    void setRect(const GRect& r) override {}
    GColor onGetColor() override { return fColor; }
    void onSetColor(const GColor& c) override { fColor = c; }

    GClick* findClick(GPoint p, GWindow* wind) override {
        if (GClick* click = Shape::findClick(p, wind)) {
            return click;
        }

        int index = -1;
        for (int i = 0; i < 4; ++i) {
            if (hit_test(p.x, p.y, fPts[i].x, fPts[i].y)) {
                index = i;
            }
        }
        if (index >= 0) {
            return new GClick(p, [this, wind, index](GClick* click) {
                fPts[index] = click->curr();
                wind->requestDraw();
            });
        }
        return nullptr;
    }
private:
    GColor  fColor;
};

