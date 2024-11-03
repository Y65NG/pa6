/*
 *  Copyright 2018 Mike Reed
 */

#include "../include/GPathBuilder.h"
#include "../include/GMatrix.h"

void GPathBuilder::reset() {
    fPts.clear();
    fVbs.clear();
}

void GPathBuilder::moveTo(GPoint p) {
    fPts.push_back(p);
    fVbs.push_back(GPathVerb::kMove);
}

void GPathBuilder::lineTo(GPoint p) {
    assert(fVbs.size() > 0);
    fPts.push_back(p);
    fVbs.push_back(GPathVerb::kLine);
}

void GPathBuilder::quadTo(GPoint p1, GPoint p2) {
    assert(fVbs.size() > 0);
    fPts.push_back(p1);
    fPts.push_back(p2);
    fVbs.push_back(GPathVerb::kQuad);
}

void GPathBuilder::cubicTo(GPoint p1, GPoint p2, GPoint p3) {
    assert(fVbs.size() > 0);
    fPts.push_back(p1);
    fPts.push_back(p2);
    fPts.push_back(p3);
    fVbs.push_back(GPathVerb::kCubic);
}

void GPathBuilder::transform(const GMatrix& m) {
    m.mapPoints(fPts.data(), fPts.size());
}

std::shared_ptr<GPath> GPathBuilder::detach() {
    auto path = std::make_shared<GPath>(std::move(fPts), std::move(fVbs));
    this->reset();
    return path;
}

/////////////////////////////////////////////////////////////

static bool is_identity(const GMatrix& m) {
    return m[0] == 1 && m[3] == 1 &&
           m[1] == 0 && m[2] == 0 && m[4] == 0 && m[5] == 0;
}

std::shared_ptr<GPath> GPath::transform(const GMatrix& m) const {
    if (fPts.empty() || is_identity(m)) {
        return const_cast<GPath*>(this)->shared_from_this();
    }
    std::vector<GPoint> dst(fPts.size());
    m.mapPoints(dst.data(), fPts.data(), fPts.size());
    return std::make_shared<GPath>(std::move(dst), fVbs);
}

GPath::Iter::Iter(const GPath& path) {
    fCurrPt = path.fPts.data();
    fCurrVb = path.fVbs.data();
    fStopVb = fCurrVb + path.fVbs.size();
}

std::optional<GPathVerb> GPath::Iter::next(GPoint pts[]) {
    assert(fCurrVb <= fStopVb);
    if (fCurrVb == fStopVb) {
        return {};
    }
    GPathVerb v = *fCurrVb++;
    switch (v) {
        case kMove:
            pts[0] = *fCurrPt++;
            break;
        case kLine:
            pts[0] = fCurrPt[-1];
            pts[1] = *fCurrPt++;
            break;
        case kQuad:
            pts[0] = fCurrPt[-1];
            pts[1] = *fCurrPt++;
            pts[2] = *fCurrPt++;
            break;
        case kCubic:
            pts[0] = fCurrPt[-1];
            pts[1] = *fCurrPt++;
            pts[2] = *fCurrPt++;
            pts[3] = *fCurrPt++;
            break;
    }
    return v;
}

constexpr int kDoneVerb = -1;

GPath::Edger::Edger(const GPath& path) {
    fPrevMove = nullptr;
    fCurrPt = path.fPts.data();
    fCurrVb = path.fVbs.data();
    fStopVb = fCurrVb + path.fVbs.size();
    fPrevVerb = kDoneVerb;
}

std::optional<GPathVerb> GPath::Edger::next(GPoint pts[]) {
    assert(fCurrVb <= fStopVb);
    bool do_return = false;
    while (fCurrVb < fStopVb) {
        switch (*fCurrVb++) {
            case kMove:
                if (fPrevVerb == kLine) {
                    pts[0] = fCurrPt[-1];
                    pts[1] = *fPrevMove;
                    do_return = true;
                }
                fPrevMove = fCurrPt++;
                fPrevVerb = kMove;
                break;
            case kLine:
                pts[0] = fCurrPt[-1];
                pts[1] = *fCurrPt++;
                fPrevVerb = kLine;
                return kLine;
            case kQuad:
                pts[0] = fCurrPt[-1];
                pts[1] = *fCurrPt++;
                pts[2] = *fCurrPt++;
                fPrevVerb = kQuad;
                return kQuad;
            case kCubic:
                pts[0] = fCurrPt[-1];
                pts[1] = *fCurrPt++;
                pts[2] = *fCurrPt++;
                pts[3] = *fCurrPt++;
                fPrevVerb = kCubic;
                return kCubic;
        }
        if (do_return) {
            return kLine;
        }
    }
    if (fPrevVerb >= kLine && fPrevVerb <= kCubic) {
        pts[0] = fCurrPt[-1];
        pts[1] = *fPrevMove;
        fPrevVerb = kDoneVerb;
        return kLine;
    } else {
        return {};
    }
}
