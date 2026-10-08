package main

import (
	"image"
	"image/color"
)

// A 5x7 bitmap font for the contact-sheet label "Page N" (the standard
// library has no text drawing). Each glyph is seven rows of five columns.
var glyphs = map[rune][7]string{
	'P': {"####.", "#...#", "#...#", "####.", "#....", "#....", "#...."},
	'a': {".....", ".....", ".###.", "....#", ".####", "#...#", ".####"},
	'g': {".....", ".####", "#...#", "#...#", ".####", "....#", ".###."},
	'e': {".....", ".....", ".###.", "#...#", "#####", "#....", ".###."},
	'0': {".###.", "#...#", "#..##", "#.#.#", "##..#", "#...#", ".###."},
	'1': {"..#..", ".##..", "..#..", "..#..", "..#..", "..#..", ".###."},
	'2': {".###.", "#...#", "....#", "...#.", "..#..", ".#...", "#####"},
	'3': {".###.", "#...#", "....#", "..##.", "....#", "#...#", ".###."},
	'4': {"...#.", "..##.", ".#.#.", "#..#.", "#####", "...#.", "...#."},
	'5': {"#####", "#....", "####.", "....#", "....#", "#...#", ".###."},
	'6': {"..##.", ".#...", "#....", "####.", "#...#", "#...#", ".###."},
	'7': {"#####", "....#", "...#.", "..#..", ".#...", ".#...", ".#..."},
	'8': {".###.", "#...#", "#...#", ".###.", "#...#", "#...#", ".###."},
	'9': {".###.", "#...#", "#...#", ".####", "....#", "...#.", ".##.."},
}

const glyphScale = 2

// drawLabel draws black text with its top-left corner at (x0, y0).
func drawLabel(dst *image.RGBA, x0, y0 int, text string) {
	x := x0
	for _, r := range text {
		g, ok := glyphs[r]
		if ok {
			for row, line := range g {
				for col, c := range line {
					if c != '#' {
						continue
					}
					for dy := 0; dy < glyphScale; dy++ {
						for dx := 0; dx < glyphScale; dx++ {
							dst.SetRGBA(x+col*glyphScale+dx, y0+row*glyphScale+dy, color.RGBA{0, 0, 0, 255})
						}
					}
				}
			}
		}
		x += 6 * glyphScale
	}
}
