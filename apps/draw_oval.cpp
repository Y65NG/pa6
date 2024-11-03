/**
 *  Copyright 2018 Mike Reed
 */

#include "../include/GPathBuilder.h"

class OvalShape : public Shape {
    GRect fRect;
    GColor  fColor;

    std::shared_ptr<GPath> makePath() const {
        GPathBuilder bu;
        bu.addCircle({0.5f, 0.5f}, 0.5f);
        bu.transform(GMatrix::Translate(fRect.left, fRect.top)
                     * GMatrix::Scale(fRect.width(), fRect.height()));
        return bu.detach();
    }

    void drawPts(GCanvas* canvas, const GPoint p[], int count) const {
        GPaint paint({0,0,0,0.5f});
        const float r = 2.5f;
        for (int i = 0; i < count; ++i) {
            canvas->drawRect({p[i].x - r, p[i].y - r, p[i].x + r, p[i].y + r}, paint);
        }
    }
public:
    OvalShape() {
        fRect = GRect::LTRB(40, 70, 200, 200);
        fColor = { 1, 0, 0, 1 };
    }

    void onDraw(GCanvas* canvas, const GPaint& paint) override {
        canvas->drawPath(*this->makePath(), paint);
    }

    void drawHilite(GCanvas* canvas) override {
        this->Shape::drawHilite(canvas);

        auto path = this->makePath();
        GPath::Iter iter(*path);
        GPoint pts[GPath::kMaxNextPoints];
        while (auto verb = iter.next(pts)) {
            switch (*verb) {
                case GPathVerb::kMove:  this->drawPts(canvas, pts, 1); break;
                case GPathVerb::kLine:  this->drawPts(canvas, pts, 1); break;
                case GPathVerb::kQuad:  this->drawPts(canvas, pts, 2); break;
                case GPathVerb::kCubic: this->drawPts(canvas, pts, 3); break;
            }
        }
    }

    GRect getRect() override { return fRect; }

    void setRect(const GRect& r) override { fRect = r; }
    GColor onGetColor() override { return fColor; }
    void onSetColor(const GColor& c) override { fColor = c; }
};

