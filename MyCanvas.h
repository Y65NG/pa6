/*
 *  Copyright 2024 <me>
 */

#ifndef _g_starter_canvas_h_
#define _g_starter_canvas_h_

#include "include/GBitmap.h"
#include "include/GCanvas.h"
#include "include/GColor.h"
#include "include/GMatrix.h"
#include "include/GPaint.h"
#include "include/GRect.h"
#include <vector>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {}

    void clear(const GColor&) override;
    void drawRect(const GRect&, const GPaint&) override;
    void drawConvexPolygon(const GPoint[], int, const GPaint&) override;
    // void fillRect(const GRect& rect, const GColor& color);
    void save() override;
    void restore() override;
    void concat(const GMatrix&) override;
    GMatrix getCtm();
    void drawPath(const GPath&, const GPaint&) override;
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count,
                  const int indices[], const GPaint&) override;
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level,
                  const GPaint&) override;

private:
    // Note: we store a copy of the bitmap
    const GBitmap fDevice;

    // Add whatever other fields you need
    GMatrix ctm = GMatrix();
    // GMatrix lastCtm = GMatrix();
    std::vector<GMatrix> copies;
    // std::vector<GMatrix> matrices;
};

#endif
