package main

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"math"
	"os"
)

func loadRGBA(path string) (*image.RGBA, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	src, err := png.Decode(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	b := src.Bounds()
	dst := image.NewRGBA(image.Rect(0, 0, b.Dx(), b.Dy()))
	draw.Draw(dst, dst.Bounds(), src, b.Min, draw.Src)
	return dst, nil
}

func savePNG(path string, im image.Image) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := png.Encode(f, im); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// thumbSize is the size Pillow's Image.thumbnail((maxW, maxH)) produces: an
// image that already fits is unchanged, otherwise the aspect ratio is kept.
func thumbSize(w, h, maxW, maxH int) (int, int) {
	if w <= maxW && h <= maxH {
		return w, h
	}
	aspect := float64(w) / float64(h)
	x, y := float64(maxW), float64(maxH)
	roundAspect := func(n float64, key func(float64) float64) int {
		lo, hi := math.Floor(n), math.Ceil(n)
		best := lo
		if key(hi) < key(lo) {
			best = hi
		}
		if best < 1 {
			best = 1
		}
		return int(best)
	}
	if x/y >= aspect {
		xi := roundAspect(y*aspect, func(n float64) float64 { return math.Abs(aspect - n/y) })
		return xi, maxH
	}
	yi := roundAspect(x/aspect, func(n float64) float64 {
		if n == 0 {
			return 0
		}
		return math.Abs(aspect - x/n)
	})
	return maxW, yi
}

func bicubic(x float64) float64 {
	const a = -0.5
	x = math.Abs(x)
	switch {
	case x < 1:
		return ((a+2)*x-(a+3))*x*x + 1
	case x < 2:
		return (((x-5)*x+8)*x - 4) * a
	}
	return 0
}

type tap struct {
	min int
	w   []float64
}

// taps computes Pillow-style (antialiasing) convolution weights for resizing
// a line of in samples to out samples with the bicubic filter.
func taps(in, out int) []tap {
	scale := float64(in) / float64(out)
	filterScale := math.Max(scale, 1)
	support := 2 * filterScale
	res := make([]tap, out)
	for xx := 0; xx < out; xx++ {
		center := (float64(xx) + 0.5) * scale
		lo := int(center - support + 0.5)
		if lo < 0 {
			lo = 0
		}
		hi := int(center + support + 0.5)
		if hi > in {
			hi = in
		}
		w := make([]float64, hi-lo)
		sum := 0.0
		for k := range w {
			w[k] = bicubic((float64(k+lo) - center + 0.5) / filterScale)
			sum += w[k]
		}
		for k := range w {
			w[k] /= sum
		}
		res[xx] = tap{lo, w}
	}
	return res
}

func clip8(v float64) uint8 {
	v = math.Floor(v + 0.5)
	if v < 0 {
		return 0
	}
	if v > 255 {
		return 255
	}
	return uint8(v)
}

// resize scales src to w x h: horizontal pass then vertical pass, 8-bit
// rounding in between, as Pillow does for RGB images.
func resize(src *image.RGBA, w, h int) *image.RGBA {
	sw, sh := src.Bounds().Dx(), src.Bounds().Dy()
	if w == sw && h == sh {
		return src
	}
	tmp := image.NewRGBA(image.Rect(0, 0, w, sh))
	hx := taps(sw, w)
	for y := 0; y < sh; y++ {
		for x := 0; x < w; x++ {
			var acc [3]float64
			for k, wt := range hx[x].w {
				o := src.PixOffset(hx[x].min+k, y)
				for c := 0; c < 3; c++ {
					acc[c] += wt * float64(src.Pix[o+c])
				}
			}
			o := tmp.PixOffset(x, y)
			tmp.Pix[o], tmp.Pix[o+1], tmp.Pix[o+2], tmp.Pix[o+3] = clip8(acc[0]), clip8(acc[1]), clip8(acc[2]), 255
		}
	}
	dst := image.NewRGBA(image.Rect(0, 0, w, h))
	vy := taps(sh, h)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			var acc [3]float64
			for k, wt := range vy[y].w {
				o := tmp.PixOffset(x, vy[y].min+k)
				for c := 0; c < 3; c++ {
					acc[c] += wt * float64(tmp.Pix[o+c])
				}
			}
			o := dst.PixOffset(x, y)
			dst.Pix[o], dst.Pix[o+1], dst.Pix[o+2], dst.Pix[o+3] = clip8(acc[0]), clip8(acc[1]), clip8(acc[2]), 255
		}
	}
	return dst
}

var tileGrey = color.RGBA{0xdd, 0xdd, 0xdd, 255}

// makeTile is a 430x560 grey tile with the page thumbnail (fitted into
// 408x530) centred horizontally 22 px from the top and the label "Page N" at
// (8, 4).
func makeTile(page *image.RGBA, number int) *image.RGBA {
	w, h := thumbSize(page.Bounds().Dx(), page.Bounds().Dy(), 408, 530)
	thumb := resize(page, w, h)
	tile := image.NewRGBA(image.Rect(0, 0, tileW, tileH))
	draw.Draw(tile, tile.Bounds(), image.NewUniform(tileGrey), image.Point{}, draw.Src)
	at := image.Pt((tileW-w)/2, 22)
	draw.Draw(tile, image.Rectangle{Min: at, Max: at.Add(image.Pt(w, h))}, thumb, image.Point{}, draw.Src)
	drawLabel(tile, 8, 4, fmt.Sprintf("Page %d", number))
	return tile
}

// makeSheet lays tiles out three per row on a white background.
func makeSheet(tiles []*image.RGBA) *image.RGBA {
	rows := (len(tiles) + tilesPerRow - 1) / tilesPerRow
	sheet := image.NewRGBA(image.Rect(0, 0, tilesPerRow*tileW, tileH*rows))
	draw.Draw(sheet, sheet.Bounds(), image.NewUniform(color.White), image.Point{}, draw.Src)
	for i, t := range tiles {
		at := image.Pt((i%tilesPerRow)*tileW, (i/tilesPerRow)*tileH)
		draw.Draw(sheet, image.Rectangle{Min: at, Max: at.Add(image.Pt(tileW, tileH))}, t, image.Point{}, draw.Src)
	}
	return sheet
}
