Clazz.declarePackage ("J.g3d");
Clazz.load (["J.api.JmolRendererInterface", "J.util.GData", "JU.V3", "J.util.Normix"], "J.g3d.Graphics3D", ["java.lang.NullPointerException", "java.util.Arrays", "javajs.awt.Font", "JU.AU", "J.api.Interface", "J.constant.EnumStereoMode", "J.g3d.CylinderRenderer", "$.ImageRenderer", "$.LineRenderer", "$.Pixelator", "$.PixelatorShaded", "$.Platform3D", "$.SphereRenderer", "$.TextRenderer", "$.TextSorter", "$.TextString", "J.util.C", "$.Shader"], function () {
c$ = Clazz.decorateAsClass (function () {
this.platform = null;
this.line3d = null;
this.sphere3d = null;
this.cylinder3d = null;
this.triangle3d = null;
this.circle3d = null;
this.hermite3d = null;
this.isFullSceneAntialiasingEnabled = false;
this.antialias2 = false;
this.strings = null;
this.stringCount = 0;
this.anaglyphChannelBytes = null;
this.twoPass = false;
this.addAllPixels = false;
this.$haveTranslucentObjects = false;
this.pbuf = null;
this.pbufT = null;
this.zbuf = null;
this.zbufT = null;
this.translucencyMask = 0;
this.renderLow = false;
this.shadesCurrent = null;
this.anaglyphLength = 0;
this.isScreened = false;
this.argbNoisyUp = 0;
this.argbNoisyDn = 0;
this.currentFont = null;
this.pixel = null;
this.zMargin = 0;
this.aobuf = null;
this.currentShadeIndex = 0;
this.lastRawColor = 0;
this.translucencyLog = 0;
this.saveAmbient = 0;
this.saveDiffuse = 0;
this.$currentlyRendering = false;
this.vectorAB = null;
this.vectorAC = null;
this.vectorNormal = null;
this.transformedVectors = null;
this.shadeIndexes = null;
this.shadeIndexes2Sided = null;
Clazz.instantialize (this, arguments);
}, J.g3d, "Graphics3D", J.util.GData, J.api.JmolRendererInterface);
Clazz.prepareFields (c$, function () {
this.vectorAB =  new JU.V3 ();
this.vectorAC =  new JU.V3 ();
this.vectorNormal =  new JU.V3 ();
this.transformedVectors =  new Array (J.g3d.Graphics3D.normixCount);
{
for (var i = J.g3d.Graphics3D.normixCount; --i >= 0; ) this.transformedVectors[i] =  new JU.V3 ();

}this.shadeIndexes =  Clazz.newByteArray (J.g3d.Graphics3D.normixCount, 0);
this.shadeIndexes2Sided =  Clazz.newByteArray (J.g3d.Graphics3D.normixCount, 0);
});
$_V(c$, "clear", 
function () {
this.stringCount = 0;
this.strings = null;
J.g3d.TextRenderer.clearFontCache ();
});
$_V(c$, "destroy", 
function () {
this.releaseBuffers ();
this.platform = null;
this.graphicsForMetrics = null;
});
$_V(c$, "getGData", 
function () {
return this;
});
$_M(c$, "setZMargin", 
function (dz) {
this.zMargin = dz;
}, "~N");
Clazz.makeConstructor (c$, 
function () {
Clazz.superConstructor (this, J.g3d.Graphics3D, []);
});
$_V(c$, "initialize", 
function (apiPlatform) {
this.apiPlatform = apiPlatform;
this.platform =  new J.g3d.Platform3D (apiPlatform);
this.graphicsForMetrics = this.platform.getGraphicsForMetrics ();
this.line3d =  new J.g3d.LineRenderer (this);
this.sphere3d =  new J.g3d.SphereRenderer (this);
this.cylinder3d =  new J.g3d.CylinderRenderer (this);
}, "javajs.api.GenericPlatform");
$_V(c$, "addRenderer", 
function (tok) {
switch (tok) {
case 1073741880:
if (this.circle3d == null) this.circle3d = this.getRenderer ("Circle");
break;
case 553648147:
if (this.hermite3d == null) this.hermite3d = this.getRenderer ("Hermite");
case 1073742182:
if (this.triangle3d == null) this.triangle3d = this.getRenderer ("Triangle");
break;
}
}, "~N");
$_M(c$, "getRenderer", 
($fz = function (type) {
return (J.api.Interface.getOptionInterface ("g3d." + type + "Renderer")).set (this);
}, $fz.isPrivate = true, $fz), "~S");
$_V(c$, "currentlyRendering", 
function () {
return this.$currentlyRendering;
});
$_V(c$, "setWindowParameters", 
function (width, height, antialias) {
this.setWinParams (width, height, antialias);
if (this.$currentlyRendering) this.endRendering ();
}, "~N,~N,~B");
$_V(c$, "checkTranslucent", 
function (isAlphaTranslucent) {
if (isAlphaTranslucent) this.$haveTranslucentObjects = true;
return (!this.twoPass || this.twoPass && (this.$isPass2 == isAlphaTranslucent));
}, "~B");
$_V(c$, "beginRendering", 
function (rotationMatrix, translucentMode, isImageWrite, renderLow) {
if (this.$currentlyRendering) this.endRendering ();
this.renderLow = renderLow;
if (this.windowWidth != this.newWindowWidth || this.windowHeight != this.newWindowHeight || this.newAntialiasing != this.isFullSceneAntialiasingEnabled) {
this.windowWidth = this.newWindowWidth;
this.windowHeight = this.newWindowHeight;
this.isFullSceneAntialiasingEnabled = this.newAntialiasing;
this.releaseBuffers ();
}this.setRotationMatrix (rotationMatrix);
this.antialiasEnabled = this.antialiasThisFrame = this.newAntialiasing;
this.$currentlyRendering = true;
if (this.strings != null) for (var i = Math.min (this.strings.length, this.stringCount); --i >= 0; ) this.strings[i] = null;

this.stringCount = 0;
this.twoPass = true;
this.$isPass2 = false;
this.colixCurrent = 0;
this.$haveTranslucentObjects = false;
this.translucentCoverOnly = !translucentMode;
this.addAllPixels = true;
if (this.pbuf == null) {
this.platform.allocateBuffers (this.windowWidth, this.windowHeight, this.antialiasThisFrame, isImageWrite);
this.pbuf = this.platform.pBuffer;
this.zbuf = this.platform.zBuffer;
this.aobuf = null;
}this.setWidthHeight (this.antialiasThisFrame);
this.platform.clearBuffer ();
if (this.backgroundImage != null) this.plotImage (-2147483648, 0, -2147483648, this.backgroundImage, null, 0, 0, 0);
}, "JU.M3,~B,~B,~B");
$_V(c$, "setBackgroundTransparent", 
function (TF) {
if (this.platform != null) this.platform.setBackgroundTransparent (TF);
}, "~B");
$_M(c$, "releaseBuffers", 
($fz = function () {
this.pbuf = null;
this.zbuf = null;
this.pbufT = null;
this.zbufT = null;
this.aobuf = null;
this.platform.releaseBuffers ();
this.line3d.clearLineCache ();
}, $fz.isPrivate = true, $fz));
$_V(c$, "setPass2", 
function (antialiasTranslucent) {
if (!this.$haveTranslucentObjects || !this.$currentlyRendering) return false;
this.$isPass2 = true;
this.colixCurrent = 0;
this.addAllPixels = true;
if (this.pbufT == null || this.antialias2 != antialiasTranslucent) {
this.platform.allocateTBuffers (antialiasTranslucent);
this.pbufT = this.platform.pBufferT;
this.zbufT = this.platform.zBufferT;
}this.antialias2 = antialiasTranslucent;
if (this.antialiasThisFrame && !this.antialias2) this.downsampleFullSceneAntialiasing (true);
this.platform.clearTBuffer ();
return true;
}, "~B");
$_V(c$, "endRendering", 
function () {
if (!this.$currentlyRendering) return;
if (this.pbuf != null) {
if (this.$isPass2 && this.pbufT != null) for (var offset = this.pbufT.length; --offset >= 0; ) J.g3d.Graphics3D.mergeBufferPixel (this.pbuf, offset, this.pbufT[offset], this.bgcolor);

if (this.ambientOcclusion != 0) {
if (this.aobuf == null) this.aobuf =  Clazz.newIntArray (this.pbuf.length, 0);
 else for (var offset = this.pbuf.length; --offset >= 0; ) this.aobuf[offset] = 0;

this.shader.occludePixels (this.pbuf, this.zbuf, this.aobuf, this.width, this.height, this.ambientOcclusion);
}if (this.antialiasThisFrame) this.downsampleFullSceneAntialiasing (false);
}this.platform.setBackgroundColor (this.bgcolor);
this.platform.notifyEndOfRendering ();
this.$currentlyRendering = false;
});
c$.mergeBufferPixel = $_M(c$, "mergeBufferPixel", 
function (pbuf, offset, argbB, bgcolor) {
if (argbB == 0) return;
var argbA = pbuf[offset];
if (argbA == argbB) return;
if (argbA == 0) argbA = bgcolor;
var rbA = (argbA & 0x00FF00FF);
var gA = (argbA & 0x0000FF00);
var rbB = (argbB & 0x00FF00FF);
var gB = (argbB & 0x0000FF00);
var logAlpha = (argbB >> 24) & 0xF;
switch (logAlpha) {
case 0:
rbA = rbB;
gA = gB;
break;
case 1:
rbA = (((rbB << 2) + (rbB << 1) + rbB + rbA) >> 3) & 0x00FF00FF;
gA = (((gB << 2) + +(gB << 1) + gB + gA) >> 3) & 0x0000FF00;
break;
case 2:
rbA = (((rbB << 1) + rbB + rbA) >> 2) & 0x00FF00FF;
gA = (((gB << 1) + gB + gA) >> 2) & 0x0000FF00;
break;
case 3:
rbA = (((rbB << 2) + rbB + (rbA << 1) + rbA) >> 3) & 0x00FF00FF;
gA = (((gB << 2) + gB + (gA << 1) + gA) >> 3) & 0x0000FF00;
break;
case 4:
rbA = ((rbA + rbB) >> 1) & 0x00FF00FF;
gA = ((gA + gB) >> 1) & 0x0000FF00;
break;
case 5:
rbA = (((rbB << 1) + rbB + (rbA << 2) + rbA) >> 3) & 0x00FF00FF;
gA = (((gB << 1) + gB + (gA << 2) + gA) >> 3) & 0x0000FF00;
break;
case 6:
rbA = (((rbA << 1) + rbA + rbB) >> 2) & 0x00FF00FF;
gA = (((gA << 1) + gA + gB) >> 2) & 0x0000FF00;
break;
case 7:
rbA = (((rbA << 2) + (rbA << 1) + rbA + rbB) >> 3) & 0x00FF00FF;
gA = (((gA << 2) + (gA << 1) + gA + gB) >> 3) & 0x0000FF00;
break;
}
pbuf[offset] = 0xFF000000 | rbA | gA;
}, "~A,~N,~N,~N");
$_V(c$, "getScreenImage", 
function (isImageWrite) {
{
var obj = this.platform.bufferedImage;
if (isImageWrite) { this.releaseBuffers(); }
return obj;
}}, "~B");
$_V(c$, "applyAnaglygh", 
function (stereoMode, stereoColors) {
switch (stereoMode) {
case J.constant.EnumStereoMode.REDCYAN:
this.applyCyanAnaglyph ();
break;
case J.constant.EnumStereoMode.CUSTOM:
this.applyCustomAnaglyph (stereoColors);
break;
case J.constant.EnumStereoMode.REDBLUE:
this.applyBlueAnaglyph ();
break;
case J.constant.EnumStereoMode.REDGREEN:
this.applyGreenAnaglyph ();
break;
case J.constant.EnumStereoMode.DOUBLE:
break;
case J.constant.EnumStereoMode.NONE:
break;
}
}, "J.constant.EnumStereoMode,~A");
$_V(c$, "snapshotAnaglyphChannelBytes", 
function () {
if (this.$currentlyRendering) throw  new NullPointerException ();
this.anaglyphLength = this.windowWidth * this.windowHeight;
if (this.anaglyphChannelBytes == null || this.anaglyphChannelBytes.length != this.anaglyphLength) this.anaglyphChannelBytes =  Clazz.newByteArray (this.anaglyphLength, 0);
for (var i = this.anaglyphLength; --i >= 0; ) this.anaglyphChannelBytes[i] = this.pbuf[i];

});
$_M(c$, "applyCustomAnaglyph", 
function (stereoColors) {
var color1 = stereoColors[0];
var color2 = stereoColors[1] & 0x00FFFFFF;
for (var i = this.anaglyphLength; --i >= 0; ) {
var a = this.anaglyphChannelBytes[i] & 0x000000FF;
a = (a | ((a | (a << 8)) << 8)) & color2;
this.pbuf[i] = (this.pbuf[i] & color1) | a;
}
}, "~A");
$_M(c$, "applyGreenAnaglyph", 
function () {
for (var i = this.anaglyphLength; --i >= 0; ) {
var green = (this.anaglyphChannelBytes[i] & 0x000000FF) << 8;
this.pbuf[i] = (this.pbuf[i] & 0xFFFF0000) | green;
}
});
$_M(c$, "applyBlueAnaglyph", 
function () {
for (var i = this.anaglyphLength; --i >= 0; ) {
var blue = this.anaglyphChannelBytes[i] & 0x000000FF;
this.pbuf[i] = (this.pbuf[i] & 0xFFFF0000) | blue;
}
});
$_M(c$, "applyCyanAnaglyph", 
function () {
for (var i = this.anaglyphLength; --i >= 0; ) {
var blue = this.anaglyphChannelBytes[i] & 0x000000FF;
var cyan = (blue << 8) | blue;
this.pbuf[i] = this.pbuf[i] & 0xFFFF0000 | cyan;
}
});
$_V(c$, "releaseScreenImage", 
function () {
this.platform.clearScreenBufferThreaded ();
});
$_V(c$, "haveTranslucentObjects", 
function () {
return this.$haveTranslucentObjects;
});
$_M(c$, "setTempZSlab", 
function (zSlab) {
this.zSlab = zSlab;
}, "~N");
$_V(c$, "setZShade", 
function (zShade, zSlab, zDepth, zShadePower) {
if (zShade) {
this.setZShade2 (zSlab, zDepth, zShadePower);
this.pixel =  new J.g3d.PixelatorShaded (this);
} else {
this.pixel =  new J.g3d.Pixelator (this);
}}, "~B,~N,~N,~N");
$_M(c$, "downsampleFullSceneAntialiasing", 
($fz = function (downsampleZBuffer) {
var width4 = this.width;
var offset1 = 0;
var offset4 = 0;
var bgcheck = this.bgcolor;
if (downsampleZBuffer) bgcheck += ((bgcheck & 0xFF) == 0xFF ? -1 : 1);
for (var i = this.pbuf.length; --i >= 0; ) if (this.pbuf[i] == 0) this.pbuf[i] = bgcheck;

bgcheck &= 0xFFFFFF;
for (var i = this.windowHeight; --i >= 0; offset4 += width4) for (var j = this.windowWidth; --j >= 0; ++offset1) {
var argb = ((this.pbuf[offset4] >> 2) & 0x3F3F3F3F) + ((this.pbuf[offset4++ + width4] >> 2) & 0x3F3F3F3F) + ((this.pbuf[offset4] >> 2) & 0x3F3F3F3F) + ((this.pbuf[offset4++ + width4] >> 2) & 0x3F3F3F3F);
argb += (argb & 0xC0C0C0C0) >> 6;
{
this.pbuf[offset1] = argb & 0x00FFFFFF | 0xFF000000;
}}

if (downsampleZBuffer) {
offset1 = offset4 = 0;
for (var i = this.windowHeight; --i >= 0; offset4 += width4) for (var j = this.windowWidth; --j >= 0; ++offset1, ++offset4) {
var z = Math.min (this.zbuf[offset4], this.zbuf[offset4 + width4]);
z = Math.min (z, this.zbuf[++offset4]);
z = Math.min (z, this.zbuf[offset4 + width4]);
if (z != 2147483647) z >>= 1;
this.zbuf[offset1] = (this.pbuf[offset1] == bgcheck ? 2147483647 : z);
}

this.antialiasThisFrame = false;
this.setWidthHeight (false);
}}, $fz.isPrivate = true, $fz), "~B");
$_M(c$, "hasContent", 
function () {
return this.platform.hasContent ();
});
$_V(c$, "setColor", 
function (argb) {
this.argbCurrent = this.argbNoisyUp = this.argbNoisyDn = argb;
}, "~N");
$_V(c$, "setColix", 
function (colix) {
var isLast = J.util.C.isColixLastAvailable (colix);
if (!isLast && colix == this.colixCurrent && this.currentShadeIndex == -1) return true;
var mask = colix & 30720;
if (mask == 16384) return false;
if (this.renderLow) mask = 0;
var isTranslucent = (mask != 0);
this.isScreened = isTranslucent && mask == 30720;
if (!this.checkTranslucent (isTranslucent && !this.isScreened)) return false;
this.addAllPixels = this.$isPass2 || !isTranslucent;
if (this.$isPass2) {
this.translucencyMask = (mask << 13) | 0xFFFFFF;
this.translucencyLog = mask >> 11;
} else {
this.translucencyLog = 0;
}this.colixCurrent = colix;
if (isLast) {
if (this.argbCurrent != this.lastRawColor) {
if (this.argbCurrent == 0) this.argbCurrent = 0xFFFFFFFF;
this.lastRawColor = this.argbCurrent;
this.shader.setLastColix (this.argbCurrent, this.inGreyscaleMode);
}}this.shadesCurrent = this.getShades (colix);
this.currentShadeIndex = -1;
this.setColor (this.getColorArgbOrGray (colix));
return true;
}, "~N");
$_M(c$, "addPixel", 
function (offset, z, p) {
this.pixel.addPixel (offset, z, p);
}, "~N,~N,~N");
$_M(c$, "clearPixel", 
function (offset, z) {
this.pixel.clearPixel (offset, z);
}, "~N,~N");
$_V(c$, "drawFilledCircle", 
function (colixRing, colixFill, diameter, x, y, z) {
if (this.isClippedZ (z)) return;
var r = Clazz.doubleToInt ((diameter + 1) / 2);
var isClipped = x < r || x + r >= this.width || y < r || y + r >= this.height;
if (isClipped && this.isClippedXY (diameter, x, y)) return;
if (colixRing != 0 && this.setColix (colixRing)) {
if (isClipped) (this.circle3d).plotCircleCenteredClipped (x, y, z, diameter);
 else (this.circle3d).plotCircleCenteredUnclipped (x, y, z, diameter);
}if (colixFill != 0 && this.setColix (colixFill)) {
if (isClipped) (this.circle3d).plotFilledCircleCenteredClipped (x, y, z, diameter);
 else (this.circle3d).plotFilledCircleCenteredUnclipped (x, y, z, diameter);
}}, "~N,~N,~N,~N,~N,~N");
$_V(c$, "volumeRender4", 
function (diameter, x, y, z) {
if (diameter == 1) {
this.plotPixelClippedXYZ (x, y, z);
return;
}if (this.isClippedZ (z)) return;
var r = Clazz.doubleToInt ((diameter + 1) / 2);
var isClipped = x < r || x + r >= this.width || y < r || y + r >= this.height;
if (isClipped && this.isClippedXY (diameter, x, y)) return;
if (isClipped) (this.circle3d).plotFilledCircleCenteredClipped (x, y, z, diameter);
 else (this.circle3d).plotFilledCircleCenteredUnclipped (x, y, z, diameter);
}, "~N,~N,~N,~N");
$_V(c$, "fillSphereXYZ", 
function (diameter, x, y, z) {
switch (diameter) {
case 1:
this.plotPixelClippedArgb (this.argbCurrent, x, y, z);
return;
case 0:
return;
}
if (diameter <= (this.antialiasThisFrame ? 2000 : 1000)) this.sphere3d.render (this.shadesCurrent, !this.addAllPixels, diameter, x, y, z, null, null, null, -1, null, this.addAllPixels);
}, "~N,~N,~N,~N");
$_V(c$, "volumeRender", 
function (TF) {
if (TF) {
this.saveAmbient = this.shader.ambientPercent;
this.saveDiffuse = this.shader.diffusePercent;
this.setAmbientPercent (100);
this.setDiffusePercent (0);
} else {
this.setAmbientPercent (this.saveAmbient);
this.setDiffusePercent (this.saveDiffuse);
}}, "~B");
$_V(c$, "fillSphereI", 
function (diameter, center) {
this.fillSphereXYZ (diameter, center.x, center.y, center.z);
}, "~N,JU.P3i");
$_V(c$, "fillSphere", 
function (diameter, center) {
this.fillSphereXYZ (diameter, Math.round (center.x), Math.round (center.y), Math.round (center.z));
}, "~N,JU.P3");
$_V(c$, "fillEllipsoid", 
function (center, points, x, y, z, diameter, mToEllipsoidal, coef, mDeriv, selectedOctant, octantPoints) {
switch (diameter) {
case 1:
this.plotPixelClippedArgb (this.argbCurrent, x, y, z);
return;
case 0:
return;
}
if (diameter <= (this.antialiasThisFrame ? 2000 : 1000)) this.sphere3d.render (this.shadesCurrent, !this.addAllPixels, diameter, x, y, z, mToEllipsoidal, coef, mDeriv, selectedOctant, octantPoints, this.addAllPixels);
}, "JU.P3,~A,~N,~N,~N,~N,JU.M3,~A,JU.M4,~N,~A");
$_V(c$, "drawRect", 
function (x, y, z, zSlab, rWidth, rHeight) {
if (zSlab != 0 && this.isClippedZ (zSlab)) return;
var w = rWidth - 1;
var h = rHeight - 1;
var xRight = x + w;
var yBottom = y + h;
if (y >= 0 && y < this.height) this.drawHLine (x, y, z, w);
if (yBottom >= 0 && yBottom < this.height) this.drawHLine (x, yBottom, z, w);
if (x >= 0 && x < this.width) this.drawVLine (x, y, z, h);
if (xRight >= 0 && xRight < this.width) this.drawVLine (xRight, y, z, h);
}, "~N,~N,~N,~N,~N,~N");
$_M(c$, "drawHLine", 
($fz = function (x, y, z, w) {
if (w < 0) {
x += w;
w = -w;
}if (x < 0) {
w += x;
x = 0;
}if (x + w >= this.width) w = this.width - 1 - x;
var offset = x + this.width * y;
if (this.addAllPixels) {
for (var i = 0; i <= w; i++) {
if (z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset++;
}
return;
}var flipflop = ((x ^ y) & 1) != 0;
for (var i = 0; i <= w; i++) {
if ((flipflop = !flipflop) && z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset++;
}
}, $fz.isPrivate = true, $fz), "~N,~N,~N,~N");
$_M(c$, "drawVLine", 
($fz = function (x, y, z, h) {
if (h < 0) {
y += h;
h = -h;
}if (y < 0) {
h += y;
y = 0;
}if (y + h >= this.height) {
h = this.height - 1 - y;
}var offset = x + this.width * y;
if (this.addAllPixels) {
for (var i = 0; i <= h; i++) {
if (z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset += this.width;
}
return;
}var flipflop = ((x ^ y) & 1) != 0;
for (var i = 0; i <= h; i++) {
if ((flipflop = !flipflop) && z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset += this.width;
}
}, $fz.isPrivate = true, $fz), "~N,~N,~N,~N");
$_V(c$, "fillRect", 
function (x, y, z, zSlab, widthFill, heightFill) {
if (this.isClippedZ (zSlab)) return;
if (x < 0) {
widthFill += x;
if (widthFill <= 0) return;
x = 0;
}if (x + widthFill > this.width) {
widthFill = this.width - x;
if (widthFill <= 0) return;
}if (y < 0) {
heightFill += y;
if (heightFill <= 0) return;
y = 0;
}if (y + heightFill > this.height) heightFill = this.height - y;
while (--heightFill >= 0) this.plotPixelsUnclippedCount (widthFill, x, y++, z);

}, "~N,~N,~N,~N,~N,~N");
$_V(c$, "drawString", 
function (str, font3d, xBaseline, yBaseline, z, zSlab, bgColix) {
this.currentShadeIndex = 0;
if (str == null) return;
if (this.isClippedZ (zSlab)) return;
this.drawStringNoSlab (str, font3d, xBaseline, yBaseline, z, bgColix);
}, "~S,javajs.awt.Font,~N,~N,~N,~N,~N");
$_V(c$, "drawStringNoSlab", 
function (str, font3d, xBaseline, yBaseline, z, bgColix) {
if (str == null) return;
if (this.strings == null) this.strings =  new Array (10);
if (this.stringCount == this.strings.length) this.strings = JU.AU.doubleLength (this.strings);
var t =  new J.g3d.TextString ();
t.setText (str, font3d == null ? this.currentFont : (this.currentFont = font3d), this.argbCurrent, J.util.C.isColixTranslucent (bgColix) ? (this.getColorArgbOrGray (bgColix) & 0xFFFFFF) | ((bgColix & 30720) << 13) : 0, xBaseline, yBaseline, z);
this.strings[this.stringCount++] = t;
}, "~S,javajs.awt.Font,~N,~N,~N,~N");
$_V(c$, "renderAllStrings", 
function (jmolRenderer) {
if (this.strings == null) return;
if (this.stringCount >= 2) {
if (J.g3d.Graphics3D.sort == null) J.g3d.Graphics3D.sort =  new J.g3d.TextSorter ();
java.util.Arrays.sort (this.strings, J.g3d.Graphics3D.sort);
}for (var i = 0; i < this.stringCount; i++) {
var ts = this.strings[i];
this.plotText (ts.x, ts.y, ts.z, ts.argb, ts.bgargb, ts.text, ts.font, jmolRenderer);
}
this.strings = null;
this.stringCount = 0;
}, "~O");
$_V(c$, "plotText", 
function (x, y, z, argb, bgargb, text, font3d, jmolRenderer) {
J.g3d.TextRenderer.plot (x, y, z, argb, bgargb, text, font3d, this, jmolRenderer, this.antialiasThisFrame);
}, "~N,~N,~N,~N,~N,~S,javajs.awt.Font,J.api.JmolRendererInterface");
$_V(c$, "drawImage", 
function (objImage, x, y, z, zSlab, bgcolix, width, height) {
if (objImage == null || width == 0 || height == 0 || this.isClippedZ (zSlab)) return;
this.plotImage (x, y, z, objImage, null, bgcolix, width, height);
}, "~O,~N,~N,~N,~N,~N,~N,~N");
$_M(c$, "plotImage", 
function (x, y, z, image, jmolRenderer, bgcolix, width, height) {
this.setColix (bgcolix);
if (!this.$isPass2) this.translucencyMask = -1;
if (bgcolix == 0) this.argbCurrent = 0;
J.g3d.ImageRenderer.plotImage (x, y, z, image, this, jmolRenderer, this.antialiasThisFrame, this.argbCurrent, width, height);
}, "~N,~N,~N,~O,J.api.JmolRendererInterface,~N,~N,~N");
$_V(c$, "setFontFid", 
function (fid) {
this.currentFont = javajs.awt.Font.getFont3D (fid);
}, "~N");
$_V(c$, "setFont", 
function (font3d) {
this.currentFont = font3d;
}, "javajs.awt.Font");
$_V(c$, "getFont3DCurrent", 
function () {
return this.currentFont;
});
$_V(c$, "drawPixel", 
function (x, y, z) {
this.plotPixelClippedXYZ (x, y, z);
}, "~N,~N,~N");
$_V(c$, "drawPoints", 
function (count, coordinates, scale) {
if (scale > 1) {
var s2 = scale * scale * 0.8;
for (var i = -scale; i < scale; i++) {
for (var j = -scale; j < scale; j++) {
if (i * i + j * j > s2) continue;
this.plotPoints (count, coordinates, i, j);
this.plotPoints (count, coordinates, i, j);
}
}
} else {
this.plotPoints (count, coordinates, 0, 0);
}}, "~N,~A,~N");
$_V(c$, "drawDashedLine", 
function (run, rise, pointA, pointB) {
this.line3d.plotDashedLine (this.argbCurrent, !this.addAllPixels, run, rise, pointA.x, pointA.y, pointA.z, pointB.x, pointB.y, pointB.z, true);
}, "~N,~N,JU.P3i,JU.P3i");
$_V(c$, "drawDottedLine", 
function (pointA, pointB) {
this.line3d.plotDashedLine (this.argbCurrent, !this.addAllPixels, 2, 1, pointA.x, pointA.y, pointA.z, pointB.x, pointB.y, pointB.z, true);
}, "JU.P3i,JU.P3i");
$_V(c$, "drawLineXYZ", 
function (x1, y1, z1, x2, y2, z2) {
this.line3d.plotLine (this.argbCurrent, !this.addAllPixels, this.argbCurrent, !this.addAllPixels, x1, y1, z1, x2, y2, z2, true);
}, "~N,~N,~N,~N,~N,~N");
$_V(c$, "drawLine", 
function (colixA, colixB, x1, y1, z1, x2, y2, z2) {
if (!this.setColix (colixA)) colixA = 0;
var isScreenedA = !this.addAllPixels;
var argbA = this.argbCurrent;
if (!this.setColix (colixB)) colixB = 0;
if (colixA == 0 && colixB == 0) return;
this.line3d.plotLine (argbA, isScreenedA, this.argbCurrent, !this.addAllPixels, x1, y1, z1, x2, y2, z2, true);
}, "~N,~N,~N,~N,~N,~N,~N,~N");
$_V(c$, "drawLineAB", 
function (pointA, pointB) {
this.line3d.plotLine (this.argbCurrent, !this.addAllPixels, this.argbCurrent, !this.addAllPixels, pointA.x, pointA.y, pointA.z, pointB.x, pointB.y, pointB.z, true);
}, "JU.P3i,JU.P3i");
$_V(c$, "fillCylinderXYZ", 
function (colixA, colixB, endcaps, diameter, xA, yA, zA, xB, yB, zB) {
if (!this.setColix (colixA)) colixA = 0;
var isScreenedA = !this.addAllPixels;
if (!this.setColix (colixB)) colixB = 0;
if (colixA == 0 && colixB == 0) return;
this.cylinder3d.render (colixA, colixB, isScreenedA, !this.addAllPixels, endcaps, diameter, xA, yA, zA, xB, yB, zB);
}, "~N,~N,~N,~N,~N,~N,~N,~N,~N,~N");
$_V(c$, "fillCylinderScreen", 
function (endcaps, diameter, xA, yA, zA, xB, yB, zB) {
this.cylinder3d.render (this.colixCurrent, this.colixCurrent, !this.addAllPixels, !this.addAllPixels, endcaps, diameter, xA, yA, zA, xB, yB, zB);
}, "~N,~N,~N,~N,~N,~N,~N,~N");
$_V(c$, "fillCylinderScreen3I", 
function (endcaps, diameter, screenA, screenB, pt0f, pt1f, radius) {
this.cylinder3d.render (this.colixCurrent, this.colixCurrent, !this.addAllPixels, !this.addAllPixels, endcaps, diameter, screenA.x, screenA.y, screenA.z, screenB.x, screenB.y, screenB.z);
}, "~N,~N,JU.P3i,JU.P3i,JU.P3,JU.P3,~N");
$_V(c$, "fillCylinder", 
function (endcaps, diameter, screenA, screenB) {
this.cylinder3d.render (this.colixCurrent, this.colixCurrent, !this.addAllPixels, !this.addAllPixels, endcaps, diameter, screenA.x, screenA.y, screenA.z, screenB.x, screenB.y, screenB.z);
}, "~N,~N,JU.P3i,JU.P3i");
$_V(c$, "fillCylinderBits", 
function (endcaps, diameter, screenA, screenB) {
this.cylinder3d.renderBits (this.colixCurrent, this.colixCurrent, !this.addAllPixels, !this.addAllPixels, endcaps, diameter, screenA.x, screenA.y, screenA.z, screenB.x, screenB.y, screenB.z);
}, "~N,~N,JU.P3,JU.P3");
$_V(c$, "fillConeScreen", 
function (endcap, screenDiameter, screenBase, screenTip, isBarb) {
this.cylinder3d.renderCone (this.colixCurrent, !this.addAllPixels, endcap, screenDiameter, screenBase.x, screenBase.y, screenBase.z, screenTip.x, screenTip.y, screenTip.z, false, isBarb);
}, "~N,~N,JU.P3i,JU.P3i,~B");
$_V(c$, "fillConeSceen3f", 
function (endcap, screenDiameter, screenBase, screenTip) {
this.cylinder3d.renderCone (this.colixCurrent, !this.addAllPixels, endcap, screenDiameter, screenBase.x, screenBase.y, screenBase.z, screenTip.x, screenTip.y, screenTip.z, true, false);
}, "~N,~N,JU.P3,JU.P3");
$_V(c$, "drawHermite4", 
function (tension, s0, s1, s2, s3) {
(this.hermite3d).renderHermiteRope (false, tension, 0, 0, 0, s0, s1, s2, s3);
}, "~N,JU.P3i,JU.P3i,JU.P3i,JU.P3i");
$_V(c$, "drawHermite7", 
function (fill, border, tension, s0, s1, s2, s3, s4, s5, s6, s7, aspectRatio, colixBack) {
if (colixBack == 0) {
(this.hermite3d).renderHermiteRibbon (fill, border, tension, s0, s1, s2, s3, s4, s5, s6, s7, aspectRatio, 0);
return;
}(this.hermite3d).renderHermiteRibbon (fill, border, tension, s0, s1, s2, s3, s4, s5, s6, s7, aspectRatio, 1);
var colix = this.colixCurrent;
this.setColix (colixBack);
(this.hermite3d).renderHermiteRibbon (fill, border, tension, s0, s1, s2, s3, s4, s5, s6, s7, aspectRatio, -1);
this.setColix (colix);
}, "~B,~B,~N,JU.P3i,JU.P3i,JU.P3i,JU.P3i,JU.P3i,JU.P3i,JU.P3i,JU.P3i,~N,~N");
$_V(c$, "fillHermite", 
function (tension, diameterBeg, diameterMid, diameterEnd, s0, s1, s2, s3) {
(this.hermite3d).renderHermiteRope (true, tension, diameterBeg, diameterMid, diameterEnd, s0, s1, s2, s3);
}, "~N,~N,~N,~N,JU.P3i,JU.P3i,JU.P3i,JU.P3i");
$_V(c$, "drawTriangle3C", 
function (screenA, colixA, screenB, colixB, screenC, colixC, check) {
if ((check & 1) == 1) this.drawLine (colixA, colixB, screenA.x, screenA.y, screenA.z, screenB.x, screenB.y, screenB.z);
if ((check & 2) == 2) this.drawLine (colixB, colixC, screenB.x, screenB.y, screenB.z, screenC.x, screenC.y, screenC.z);
if ((check & 4) == 4) this.drawLine (colixA, colixC, screenA.x, screenA.y, screenA.z, screenC.x, screenC.y, screenC.z);
}, "JU.P3i,~N,JU.P3i,~N,JU.P3i,~N,~N");
$_V(c$, "drawTriangle3I", 
function (screenA, screenB, screenC, check) {
if ((check & 1) == 1) this.line3d.plotLine (this.argbCurrent, !this.addAllPixels, this.argbCurrent, !this.addAllPixels, screenA.x, screenA.y, screenA.z, screenB.x, screenB.y, screenB.z, true);
if ((check & 2) == 2) this.line3d.plotLine (this.argbCurrent, !this.addAllPixels, this.argbCurrent, !this.addAllPixels, screenB.x, screenB.y, screenB.z, screenC.x, screenC.y, screenC.z, true);
if ((check & 4) == 4) this.line3d.plotLine (this.argbCurrent, !this.addAllPixels, this.argbCurrent, !this.addAllPixels, screenA.x, screenA.y, screenA.z, screenC.x, screenC.y, screenC.z, true);
}, "JU.P3i,JU.P3i,JU.P3i,~N");
$_V(c$, "fillTriangleTwoSided", 
function (normix, xScreenA, yScreenA, zScreenA, xScreenB, yScreenB, zScreenB, xScreenC, yScreenC, zScreenC) {
this.setColorNoisy (this.getShadeIndex (normix));
(this.triangle3d).fillTriangleXYZ (xScreenA, yScreenA, zScreenA, xScreenB, yScreenB, zScreenB, xScreenC, yScreenC, zScreenC, false);
}, "~N,~N,~N,~N,~N,~N,~N,~N,~N,~N");
$_V(c$, "fillTriangle3f", 
function (screenA, screenB, screenC, setNoisy) {
var i = this.getShadeIndexP3 (screenA, screenB, screenC);
if (setNoisy) this.setColorNoisy (i);
 else this.setColor (this.shadesCurrent[i]);
(this.triangle3d).fillTriangleP3f (screenA, screenB, screenC, false);
}, "JU.P3,JU.P3,JU.P3,~B");
$_V(c$, "fillTriangle3i", 
function (screenA, screenB, screenC, ptA, ptB, ptC) {
(this.triangle3d).fillTriangleP3i (screenA, screenB, screenC, false);
}, "JU.P3i,JU.P3i,JU.P3i,JU.P3,JU.P3,JU.P3");
$_V(c$, "fillTriangle", 
function (screenA, colixA, normixA, screenB, colixB, normixB, screenC, colixC, normixC, factor) {
var useGouraud;
if (!this.$isPass2 && normixA == normixB && normixA == normixC && colixA == colixB && colixA == colixC) {
this.setTriangleColixAndShadeIndex (colixA, this.getShadeIndex (normixA));
useGouraud = false;
} else {
if (!this.setTriangleTranslucency (colixA, colixB, colixC)) return;
(this.triangle3d).setGouraud (this.getShades (colixA)[this.getShadeIndex (normixA)], this.getShades (colixB)[this.getShadeIndex (normixB)], this.getShades (colixC)[this.getShadeIndex (normixC)]);
useGouraud = true;
}(this.triangle3d).fillTriangleP3if (screenA, screenB, screenC, factor, useGouraud);
}, "JU.P3i,~N,~N,JU.P3i,~N,~N,JU.P3i,~N,~N,~N");
$_V(c$, "fillTriangle3CN", 
function (screenA, colixA, normixA, screenB, colixB, normixB, screenC, colixC, normixC) {
var useGouraud;
if (!this.$isPass2 && normixA == normixB && normixA == normixC && colixA == colixB && colixA == colixC) {
this.setTriangleColixAndShadeIndex (colixA, this.getShadeIndex (normixA));
useGouraud = false;
} else {
if (!this.setTriangleTranslucency (colixA, colixB, colixC)) return;
(this.triangle3d).setGouraud (this.getShades (colixA)[this.getShadeIndex (normixA)], this.getShades (colixB)[this.getShadeIndex (normixB)], this.getShades (colixC)[this.getShadeIndex (normixC)]);
useGouraud = true;
}(this.triangle3d).fillTriangleP3i (screenA, screenB, screenC, useGouraud);
}, "JU.P3i,~N,~N,JU.P3i,~N,~N,JU.P3i,~N,~N");
$_M(c$, "setTriangleColixAndShadeIndex", 
($fz = function (colix, shadeIndex) {
if (colix == this.colixCurrent && this.currentShadeIndex == shadeIndex) return;
this.currentShadeIndex = -1;
this.setColix (colix);
this.setColorNoisy (shadeIndex);
}, $fz.isPrivate = true, $fz), "~N,~N");
$_M(c$, "setTriangleTranslucency", 
($fz = function (colixA, colixB, colixC) {
if (!this.$isPass2) return true;
var maskA = colixA & 30720;
var maskB = colixB & 30720;
var maskC = colixC & 30720;
maskA &= -16385;
maskB &= -16385;
maskC &= -16385;
var mask = J.util.GData.roundInt (Clazz.doubleToInt ((maskA + maskB + maskC) / 3)) & 30720;
this.translucencyMask = (mask << 13) | 0xFFFFFF;
return true;
}, $fz.isPrivate = true, $fz), "~N,~N,~N");
$_V(c$, "drawQuadrilateral", 
function (colix, screenA, screenB, screenC, screenD) {
this.setColix (colix);
this.drawLineAB (screenA, screenB);
this.drawLineAB (screenB, screenC);
this.drawLineAB (screenC, screenD);
this.drawLineAB (screenD, screenA);
}, "~N,JU.P3i,JU.P3i,JU.P3i,JU.P3i");
$_V(c$, "fillQuadrilateral", 
function (screenA, screenB, screenC, screenD) {
this.setColorNoisy (this.getShadeIndexP3 (screenA, screenB, screenC));
(this.triangle3d).fillTriangleP3f (screenA, screenB, screenC, false);
(this.triangle3d).fillTriangleP3f (screenA, screenC, screenD, false);
}, "JU.P3,JU.P3,JU.P3,JU.P3");
$_V(c$, "fillQuadrilateral3i", 
function (screenA, colixA, normixA, screenB, colixB, normixB, screenC, colixC, normixC, screenD, colixD, normixD) {
this.fillTriangle3CN (screenA, colixA, normixA, screenB, colixB, normixB, screenC, colixC, normixC);
this.fillTriangle3CN (screenA, colixA, normixA, screenC, colixC, normixC, screenD, colixD, normixD);
}, "JU.P3i,~N,~N,JU.P3i,~N,~N,JU.P3i,~N,~N,JU.P3i,~N,~N");
$_V(c$, "drawSurface", 
function (meshSurface, colix) {
}, "J.util.MeshSurface,~N");
$_M(c$, "plotPixelClippedXYZ", 
function (x, y, z) {
if (this.isClipped3 (x, y, z)) return;
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
}, "~N,~N,~N");
$_V(c$, "plotPixelClippedP3i", 
function (screen) {
this.plotPixelClippedXYZ (screen.x, screen.y, screen.z);
}, "JU.P3i");
$_M(c$, "plotPixelClippedArgb", 
function (argb, x, y, z) {
if (this.isClipped3 (x, y, z)) return;
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.addPixel (offset, z, argb);
}, "~N,~N,~N,~N");
$_V(c$, "plotImagePixel", 
function (argb, x, y, z, shade, bgargb) {
if (this.isClipped (x, y)) return;
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.shadeTextPixel (offset, z, argb, bgargb, shade);
}, "~N,~N,~N,~N,~N,~N");
$_M(c$, "plotPixelClippedScreened", 
function (argb, isScreened, x, y, z) {
if (this.isClipped3 (x, y, z)) return;
if (isScreened && ((x ^ y) & 1) != 0) return;
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.addPixel (offset, z, argb);
}, "~N,~B,~N,~N,~N");
$_M(c$, "plotPixelUnclipped", 
function (x, y, z) {
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
}, "~N,~N,~N");
$_M(c$, "plotPixelUnclippedArgb", 
function (argb, x, y, z) {
var offset = y * this.width + x;
if (z < this.zbuf[offset]) this.addPixel (offset, z, argb);
}, "~N,~N,~N,~N");
$_M(c$, "plotPixelsClipped", 
function (count, x, y, z) {
if (y < 0 || y >= this.height || x >= this.width) return;
if (x < 0) {
count += x;
x = 0;
}if (count + x > this.width) count = this.width - x;
if (count <= 0) return;
var offsetPbuf = y * this.width + x;
var offsetMax = offsetPbuf + count;
var step = 1;
if (!this.addAllPixels) {
step = 2;
if (((x ^ y) & 1) != 0) ++offsetPbuf;
}while (offsetPbuf < offsetMax) {
if (z < this.zbuf[offsetPbuf]) this.addPixel (offsetPbuf, z, this.argbCurrent);
offsetPbuf += step;
}
}, "~N,~N,~N,~N");
$_M(c$, "plotPixelsClippedRaster", 
function (count, x, y, zAtLeft, zPastRight, rgb16Left, rgb16Right) {
if (count <= 0 || y < 0 || y >= this.height || x >= this.width || (zAtLeft < this.slab && zPastRight < this.slab) || (zAtLeft > this.depth && zPastRight > this.depth)) return;
var seed = (x << 16) + (y << 1) ^ 0x33333333;
var zScaled = (zAtLeft << 10) + (512);
var dz = zPastRight - zAtLeft;
var roundFactor = Clazz.doubleToInt (count / 2);
var zIncrementScaled = J.util.GData.roundInt (Clazz.doubleToInt (((dz << 10) + (dz >= 0 ? roundFactor : -roundFactor)) / count));
if (x < 0) {
x = -x;
zScaled += zIncrementScaled * x;
count -= x;
if (count <= 0) return;
x = 0;
}if (count + x > this.width) count = this.width - x;
var flipflop = ((x ^ y) & 1) != 0;
var offsetPbuf = y * this.width + x;
if (rgb16Left == null) {
while (--count >= 0) {
if (this.addAllPixels || (flipflop = !flipflop) == true) {
var z = zScaled >> 10;
if (z >= this.slab && z <= this.depth && z < this.zbuf[offsetPbuf]) {
seed = ((seed << 16) + (seed << 1) + seed) & 0x7FFFFFFF;
var bits = (seed >> 16) & 0x07;
this.addPixel (offsetPbuf, z, bits == 0 ? this.argbNoisyDn : (bits == 1 ? this.argbNoisyUp : this.argbCurrent));
}}++offsetPbuf;
zScaled += zIncrementScaled;
}
} else {
var rScaled = rgb16Left.rScaled << 8;
var rIncrement = Clazz.doubleToInt (((rgb16Right.rScaled - rgb16Left.rScaled) << 8) / count);
var gScaled = rgb16Left.gScaled;
var gIncrement = Clazz.doubleToInt ((rgb16Right.gScaled - gScaled) / count);
var bScaled = rgb16Left.bScaled;
var bIncrement = Clazz.doubleToInt ((rgb16Right.bScaled - bScaled) / count);
while (--count >= 0) {
if (this.addAllPixels || (flipflop = !flipflop)) {
var z = zScaled >> 10;
if (z >= this.slab && z <= this.depth && z < this.zbuf[offsetPbuf]) this.addPixel (offsetPbuf, z, 0xFF000000 | (rScaled & 0xFF0000) | (gScaled & 0xFF00) | ((bScaled >> 8) & 0xFF));
}++offsetPbuf;
zScaled += zIncrementScaled;
rScaled += rIncrement;
gScaled += gIncrement;
bScaled += bIncrement;
}
}}, "~N,~N,~N,~N,~N,J.util.Rgb16,J.util.Rgb16");
$_M(c$, "plotPixelsUnclippedRaster", 
function (count, x, y, zAtLeft, zPastRight, rgb16Left, rgb16Right) {
if (count <= 0) return;
var seed = ((x << 16) + (y << 1) ^ 0x33333333) & 0x7FFFFFFF;
var flipflop = ((x ^ y) & 1) != 0;
var zScaled = (zAtLeft << 10) + (512);
var dz = zPastRight - zAtLeft;
var roundFactor = Clazz.doubleToInt (count / 2);
var zIncrementScaled = J.util.GData.roundInt (Clazz.doubleToInt (((dz << 10) + (dz >= 0 ? roundFactor : -roundFactor)) / count));
var offsetPbuf = y * this.width + x;
if (rgb16Left == null) {
while (--count >= 0) {
if (this.addAllPixels || (flipflop = !flipflop)) {
var z = zScaled >> 10;
if (z < this.zbuf[offsetPbuf]) {
seed = ((seed << 16) + (seed << 1) + seed) & 0x7FFFFFFF;
var bits = (seed >> 16) & 0x07;
this.addPixel (offsetPbuf, z, bits == 0 ? this.argbNoisyDn : (bits == 1 ? this.argbNoisyUp : this.argbCurrent));
}}++offsetPbuf;
zScaled += zIncrementScaled;
}
} else {
var rScaled = rgb16Left.rScaled << 8;
var rIncrement = J.util.GData.roundInt (Clazz.doubleToInt (((rgb16Right.rScaled - rgb16Left.rScaled) << 8) / count));
var gScaled = rgb16Left.gScaled;
var gIncrement = J.util.GData.roundInt (Clazz.doubleToInt ((rgb16Right.gScaled - gScaled) / count));
var bScaled = rgb16Left.bScaled;
var bIncrement = J.util.GData.roundInt (Clazz.doubleToInt ((rgb16Right.bScaled - bScaled) / count));
while (--count >= 0) {
if (this.addAllPixels || (flipflop = !flipflop)) {
var z = zScaled >> 10;
if (z < this.zbuf[offsetPbuf]) this.addPixel (offsetPbuf, z, 0xFF000000 | (rScaled & 0xFF0000) | (gScaled & 0xFF00) | ((bScaled >> 8) & 0xFF));
}++offsetPbuf;
zScaled += zIncrementScaled;
rScaled += rIncrement;
gScaled += gIncrement;
bScaled += bIncrement;
}
}}, "~N,~N,~N,~N,~N,J.util.Rgb16,J.util.Rgb16");
$_M(c$, "plotPixelsUnclippedCount", 
function (count, x, y, z) {
var offsetPbuf = y * this.width + x;
if (this.addAllPixels) {
while (--count >= 0) {
if (z < this.zbuf[offsetPbuf]) this.addPixel (offsetPbuf, z, this.argbCurrent);
++offsetPbuf;
}
} else {
var offsetMax = offsetPbuf + count;
if (((x ^ y) & 1) != 0) if (++offsetPbuf == offsetMax) return;
do {
if (z < this.zbuf[offsetPbuf]) this.addPixel (offsetPbuf, z, this.argbCurrent);
offsetPbuf += 2;
} while (offsetPbuf < offsetMax);
}}, "~N,~N,~N,~N");
$_M(c$, "plotPoints", 
($fz = function (count, coordinates, xOffset, yOffset) {
for (var i = count * 3; i > 0; ) {
var z = coordinates[--i];
var y = coordinates[--i] + yOffset;
var x = coordinates[--i] + xOffset;
if (this.isClipped3 (x, y, z)) continue;
var offset = y * this.width + x++;
if (z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
if (this.antialiasThisFrame) {
offset = y * this.width + x;
if (!this.isClipped3 (x, y, z) && z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset = (++y) * this.width + x;
if (!this.isClipped3 (x, y, z) && z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
offset = y * this.width + (--x);
if (!this.isClipped3 (x, y, z) && z < this.zbuf[offset]) this.addPixel (offset, z, this.argbCurrent);
}}
}, $fz.isPrivate = true, $fz), "~N,~A,~N,~N");
$_M(c$, "setColorNoisy", 
function (shadeIndex) {
this.currentShadeIndex = shadeIndex;
this.argbCurrent = this.shadesCurrent[shadeIndex];
this.argbNoisyUp = this.shadesCurrent[shadeIndex < J.util.Shader.shadeIndexLast ? shadeIndex + 1 : J.util.Shader.shadeIndexLast];
this.argbNoisyDn = this.shadesCurrent[shadeIndex > 0 ? shadeIndex - 1 : 0];
}, "~N");
$_V(c$, "setNoisySurfaceShade", 
function (screenA, screenB, screenC) {
this.vectorAB.set (screenB.x - screenA.x, screenB.y - screenA.y, screenB.z - screenA.z);
var shadeIndex;
if (screenC == null) {
shadeIndex = this.shader.getShadeIndex (-this.vectorAB.x, -this.vectorAB.y, this.vectorAB.z);
} else {
this.vectorAC.set (screenC.x - screenA.x, screenC.y - screenA.y, screenC.z - screenA.z);
this.vectorAB.cross (this.vectorAB, this.vectorAC);
shadeIndex = this.vectorAB.z >= 0 ? this.shader.getShadeIndex (-this.vectorAB.x, -this.vectorAB.y, this.vectorAB.z) : this.shader.getShadeIndex (this.vectorAB.x, this.vectorAB.y, -this.vectorAB.z);
}if (shadeIndex > J.util.Shader.shadeIndexNoisyLimit) shadeIndex = J.util.Shader.shadeIndexNoisyLimit;
this.setColorNoisy (shadeIndex);
}, "JU.P3i,JU.P3i,JU.P3i");
$_M(c$, "getShadeIndexP3", 
($fz = function (screenA, screenB, screenC) {
this.vectorAB.sub2 (screenB, screenA);
this.vectorAC.sub2 (screenC, screenA);
this.vectorNormal.cross (this.vectorAB, this.vectorAC);
var i = (this.vectorNormal.z >= 0 ? this.shader.getShadeIndex (-this.vectorNormal.x, -this.vectorNormal.y, this.vectorNormal.z) : this.shader.getShadeIndex (this.vectorNormal.x, this.vectorNormal.y, -this.vectorNormal.z));
return i;
}, $fz.isPrivate = true, $fz), "JU.P3,JU.P3,JU.P3");
$_V(c$, "renderBackground", 
function (jmolRenderer) {
if (this.backgroundImage != null) this.plotImage (-2147483648, 0, -2147483648, this.backgroundImage, jmolRenderer, 0, 0, 0);
}, "J.api.JmolRendererInterface");
$_V(c$, "drawAtom", 
function (atom) {
this.fillSphereXYZ (atom.sD, atom.sX, atom.sY, atom.sZ);
}, "J.modelset.Atom");
$_V(c$, "getExportType", 
function () {
return 0;
});
$_V(c$, "getExportName", 
function () {
return null;
});
$_M(c$, "canDoTriangles", 
function () {
return true;
});
$_M(c$, "isCartesianExport", 
function () {
return false;
});
$_V(c$, "initializeExporter", 
function (viewer, privateKey, g3d, params) {
return null;
}, "J.viewer.Viewer,~N,J.util.GData,java.util.Map");
$_V(c$, "finalizeOutput", 
function () {
return null;
});
$_V(c$, "drawBond", 
function (atomA, atomB, colixA, colixB, endcaps, mad, bondOrder) {
}, "JU.P3,JU.P3,~N,~N,~N,~N,~N");
$_V(c$, "drawEllipse", 
function (ptAtom, ptX, ptY, fillArc, wireframeOnly) {
return false;
}, "JU.P3,JU.P3,JU.P3,~B,~B");
$_M(c$, "getPrivateKey", 
function () {
return 0;
});
$_V(c$, "clearFontCache", 
function () {
J.g3d.TextRenderer.clearFontCache ();
});
$_V(c$, "getTransformedVertexVectors", 
function () {
return this.transformedVectors;
});
$_V(c$, "isDirectedTowardsCamera", 
function (normix) {
return (normix < 0) || (this.transformedVectors[normix].z > 0);
}, "~N");
$_M(c$, "setRotationMatrix", 
function (rotationMatrix) {
var vertexVectors = J.util.Normix.getVertexVectors ();
for (var i = J.g3d.Graphics3D.normixCount; --i >= 0; ) {
var tv = this.transformedVectors[i];
rotationMatrix.rotate2 (vertexVectors[i], tv);
this.shadeIndexes[i] = this.shader.getShadeB (tv.x, -tv.y, tv.z);
this.shadeIndexes2Sided[i] = (tv.z >= 0 ? this.shadeIndexes[i] : this.shader.getShadeB (-tv.x, tv.y, -tv.z));
}
}, "JU.M3");
$_M(c$, "getShadeIndex", 
function (normix) {
return (normix == -10000 || normix == 9999 ? J.g3d.Graphics3D.nullShadeIndex : normix < 0 ? this.shadeIndexes2Sided[~normix] : this.shadeIndexes[normix]);
}, "~N");
$_V(c$, "renderCrossHairs", 
function (minMax, screenWidth, screenHeight, navOffset, navDepth) {
var antialiased = this.isAntialiased ();
this.setColix (navDepth < 0 ? 10 : navDepth > 100 ? 11 : 23);
var x = Math.max (Math.min (this.width, Math.round (navOffset.x)), 0);
var y = Math.max (Math.min (this.height, Math.round (navOffset.y)), 0);
var z = Math.round (navOffset.z) + 1;
var off = (antialiased ? 8 : 4);
var h = (antialiased ? 20 : 10);
var w = (antialiased ? 2 : 1);
this.drawRect (x - off, y, z, 0, h, w);
this.drawRect (x, y - off, z, 0, w, h);
this.drawRect (x - off, y - off, z, 0, h, h);
off = h;
h = h >> 1;
this.setColix (minMax[1] < navOffset.x ? 21 : 11);
this.drawRect (x - off, y, z, 0, h, w);
this.setColix (minMax[0] > navOffset.x ? 21 : 11);
this.drawRect (x + h, y, z, 0, h, w);
this.setColix (minMax[3] < navOffset.y ? 21 : 11);
this.drawRect (x, y - off, z, 0, w, h);
this.setColix (minMax[2] > navOffset.y ? 21 : 11);
this.drawRect (x, y + h, z, 0, w, h);
}, "~A,~N,~N,JU.P3,~N");
$_V(c$, "initializeOutput", 
function (viewer, privateKey, gdata, params) {
return false;
}, "J.viewer.Viewer,~N,J.util.GData,java.util.Map");
$_M(c$, "shadeTextPixel", 
function (offset, z, argb, bgargb, shade) {
switch (shade) {
case 8:
this.addPixel (offset, z, argb);
return;
}
if (bgargb != 0) {
J.g3d.Graphics3D.mergeBufferPixel (this.pbuf, offset, bgargb, this.bgcolor);
}shade += this.translucencyLog;
if (shade > 7) return;
J.g3d.Graphics3D.mergeBufferPixel (this.pbuf, offset, (argb & 0xFFFFFF) | shade << 24, this.bgcolor);
this.zbuf[offset] = z;
}, "~N,~N,~N,~N,~N");
Clazz.defineStatics (c$,
"sort", null);
c$.normixCount = c$.prototype.normixCount = J.util.Normix.getNormixCount ();
Clazz.defineStatics (c$,
"nullShadeIndex", 50);
});
