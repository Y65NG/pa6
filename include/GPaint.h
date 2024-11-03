/*
 *  Copyright 2016 Mike Reed
 */

#ifndef GPaint_DEFINED
#define GPaint_DEFINED

#include "GColor.h"
#include "GBlendMode.h"

class GShader;

class GPaint {
public:
    GPaint() {}
    explicit GPaint(const GColor& c) : fColor(c) {}
    explicit GPaint(std::shared_ptr<GShader> s) : fShader(std::move(s)) {}

    GColor  getColor() const { return fColor; }
    GPaint& setColor(GColor c) { fColor = c; return *this; }
    GPaint& setRGBA(float r, float g, float b, float a) {
        return this->setColor(GColor::RGBA(r, g, b, a));
    }

    float   getAlpha() const { return fColor.a; }
    GPaint& setAlpha(float alpha) {
        fColor.a = alpha;
        return *this;
    }

    GBlendMode getBlendMode() const { return fMode; }
    GPaint&    setBlendMode(GBlendMode m) { fMode = m; return *this; }

    GShader* peekShader() const { return fShader.get(); }
    std::shared_ptr<GShader> shareShader() const { return fShader; }
    GPaint&  setShader(std::shared_ptr<GShader> s) { fShader = s; return *this; }

private:
    GColor                      fColor = {0, 0, 0, 1};
    std::shared_ptr<GShader>    fShader;
    GBlendMode                  fMode = GBlendMode::kSrcOver;
};

#endif
