#ifndef GPath_DEFINED
#define GPath_DEFINED

#include "GMatrix.h"
#include "GPoint.h"
#include "GRect.h"

#include <vector>

enum GPathVerb {
    kMove,  // returns pts[0] from Iter
    kLine,  // returns pts[0]..pts[1] from Iter and Edger
    kQuad,  // returns pts[0]..pts[2] from Iter and Edger
    kCubic, // returns pts[0]..pts[3] from Iter and Edger
};

enum class GPathDirection {
    kCW,  // clockwise
    kCCW, // counter-clockwise
};

class GPath : public std::enable_shared_from_this<GPath> {
public:
    /**
     *  Return the bounds of all of the control-points in the path.
     *
     *  If there are no points, returns an empty rect (all zeros)
     */
    GRect bounds() const;

    size_t countPoints() const { return fPts.size(); }

    /**
     *  Create a new path by transforming the points in this path.
     */
    std::shared_ptr<GPath> transform(const GMatrix&) const;

    std::shared_ptr<GPath> offset(float dx, float dy) const {
        return this->transform(GMatrix::Translate(dx, dy));
    }

    // maximum number of points returned by Iter::next() and Edger::next()
    enum {
        kMaxNextPoints = 4
    };

    /**
     *  Walks the path, returning each verb that was entered, and the related point value(s)
     *      kMove : pts[0] ... the start of a new contour
     *      kLine : pts[0], pts[1] ... the start and end of that line segment
     *
     *  Make sure that the storage of pts[] has at least kMaxNextPoints entries.
     */
    class Iter {
    public:
        Iter(const GPath&);
        nonstd::optional<GPathVerb> next(GPoint pts[]);

    private:
        const GPoint*    fCurrPt;
        const GPathVerb* fCurrVb;
        const GPathVerb* fStopVb;
    };

    /**
     *  Walks the path, returning "edges" only. Thus it does not return kMove, but will return
     *  the final closing "edge" for each contour.
     *  e.g.
     *      path.moveTo(A).lineTo(B).lineTo(C).moveTo(D).lineTo(E)
     *  will return
     *      kLine   A..B
     *      kLine   B..C
     *      kLine   C..A
     *      kLine   D..E
     *      kLine   E..D
     *
     * Typical calling pattern...
     *
     *   GPoint pts[GPath::kMaxNextPoints];
     *   GPath::Edger edger(path);
     *   while (auto v = edger.next(pts)) {
     *       switch (v.value()) {
     *           case GPath::kLine: // pts[0] and pts[1] are valid
     *  }
     */
    class Edger {
    public:
        Edger(const GPath&);
        nonstd::optional<GPathVerb> next(GPoint pts[]);

    private:
        const GPoint*    fPrevMove;
        const GPoint*    fCurrPt;
        const GPathVerb* fCurrVb;
        const GPathVerb* fStopVb;
        int fPrevVerb;
    };

    /**
     *  Given 0 < t < 1, subdivide the src[] quadratic bezier at t into two new quadratics in dst[]
     *  such that
     *  0...t is stored in dst[0..2]
     *  t...1 is stored in dst[2..4]
     */
    static void ChopQuadAt(const GPoint src[3], GPoint dst[5], float t);

    /**
     *  Given 0 < t < 1, subdivide the src[] cubic bezier at t into two new cubics in dst[]
     *  such that
     *  0...t is stored in dst[0..3]
     *  t...1 is stored in dst[3..6]
     */
    static void ChopCubicAt(const GPoint src[4], GPoint dst[7], float t);

    GPath(std::vector<GPoint> pts, std::vector<GPathVerb> vbs)
        : fPts(std::move(pts))
        , fVbs(std::move(vbs))
    {}

private:
    GPath(const GPath&) = delete;
    GPath& operator=(const GPath&) = delete;

    friend class GPathBuilder;

    const std::vector<GPoint>    fPts;
    const std::vector<GPathVerb> fVbs;
};

#endif

