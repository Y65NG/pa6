#ifndef GPathBuilder_DEFINED
#define GPathBuilder_DEFINED

#include "GPoint.h"
#include "GRect.h"
#include "GPath.h"

#include <functional>
#include <optional>
#include <vector>

class GPathBuilder {
public:
    GPathBuilder() {}

    /**
     *  Erase any previously added points/verbs, restoring the path to its initial empty state.
     */
    void reset();

    /**
     *  Start a new contour at the specified coordinate.
     *  Returns a reference to this path.
     */
    void moveTo(GPoint);
    void moveTo(float x, float y) { this->moveTo({x, y}); }

    /**
     *  Connect the previous point (either from a moveTo or lineTo) with a line segment to
     *  the specified coordinate.
     *  Returns a reference to this path.
     */
    void lineTo(GPoint);
    void lineTo(float x, float y) { this->lineTo({x, y}); }

    /**
     *  Connect the previous point with a quadratic bezier to the specified coordinates.
     */
    void quadTo(GPoint, GPoint);
    void quadTo(float x1, float y1, float x2, float y2) {
        this->quadTo({x1, y1}, {x2, y2});
    }

    /**
     *  Connect the previous point with a cubic bezier to the specified coordinates.
     */
    void cubicTo(GPoint, GPoint, GPoint);
    void cubicTo(float x1, float y1, float x2, float y2, float x3, float y3) {
        this->cubicTo({x1, y1}, {x2, y2}, {x3, y3});
    }

    /**
     *  Append a new contour to this path, made up of the 4 points of the specified rect,
     *  in the specified direction.
     *
     *  In either direction the contour must begin at the top-left corner of the rect.
     */
    void addRect(const GRect&, GPathDirection = GPathDirection::kCW);

    /**
     *  Append a new contour to this path with the specified polygon.
     *  Calling this is equivalent to calling moveTo(pts[0]), lineTo(pts[1..count-1]).
     */
    void addPolygon(const GPoint pts[], int count);

    /**
     *  Append a new contour respecting the Direction. The contour should be an approximate
     *  circle (8 quadratic curves will suffice) with the specified center and radius.
     */
    void addCircle(GPoint center, float radius, GPathDirection = GPathDirection::kCW);

    void transform(const GMatrix&);

    /**
     * Return a GPath from the contents of this builder,
     * and then reset() the builder back to its empty state.
     */
    std::shared_ptr<GPath> detach();

    static std::shared_ptr<GPath> Build(std::function<void(GPathBuilder&)> proc) {
        GPathBuilder bu;
        proc(bu);
        return bu.detach();
    }

private:
    std::vector<GPoint>    fPts;
    std::vector<GPathVerb> fVbs;
};

#endif

