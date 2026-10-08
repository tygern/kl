package main

import (
	"image"
	"strings"
	"testing"
)

func TestThumbSize(t *testing.T) {
	w, h := thumbSize(510, 660, 408, 530)
	if w != 408 || h != 528 {
		t.Fatalf("got %dx%d, want 408x528", w, h)
	}
	if w, h := thumbSize(100, 100, 408, 530); w != 100 || h != 100 {
		t.Fatalf("small image resized to %dx%d", w, h)
	}
}

func TestParseInfo(t *testing.T) {
	info, err := parseInfo("Title:           A \u2013 B\nAuthor:          X\nCreationDate:    Thu Oct  8 06:24:59 2026 CDT\nPages:           12\n")
	if err != nil || info.pages != 12 || info.title == nil || *info.title != "A \u2013 B" || *info.author != "X" {
		t.Fatalf("parseInfo: %+v %v", info, err)
	}
	info, _ = parseInfo("Pages: 1\n")
	if info.title != nil {
		t.Fatal("missing title must stay nil")
	}
}

func TestParseBBox(t *testing.T) {
	text := `<doc><page width="612.000000" height="792.000000">
<word xMin="72.0" yMin="80.0" xMax="100.5" yMax="90.0">a&amp;b</word>
<word xMin="500.0" yMin="80.0" xMax="560.0" yMax="90.0">end</word>
</page><page width="612.000000" height="792.000000"></page></doc>`
	pages, err := parseBBox(text)
	if err != nil || len(pages) != 2 || len(pages[0].words) != 2 || len(pages[1].words) != 0 {
		t.Fatalf("parseBBox: %+v %v", pages, err)
	}
	if pages[0].words[0].text != "a&b" || pages[0].words[1].xMax != 560 {
		t.Fatalf("bad words: %+v", pages[0].words)
	}
}

func TestEncodeASCII(t *testing.T) {
	s := "x–y\U0001F600"
	out, err := encodeReports([]report{{PDF: "p", Title: &s, ContactSheets: []string{}}})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(out), "\"x\\u2013y\\ud83d\\ude00\"") || !strings.Contains(string(out), `"author": null`) {
		t.Fatalf("unexpected encoding:\n%s", out)
	}
	if !strings.HasSuffix(string(out), "]\n") {
		t.Fatal("missing trailing newline")
	}
}

func TestSheetLayout(t *testing.T) {
	page := image.NewRGBA(image.Rect(0, 0, 510, 660))
	tile := makeTile(page, 7)
	if tile.Bounds().Dx() != 430 || tile.Bounds().Dy() != 560 {
		t.Fatalf("tile size %v", tile.Bounds())
	}
	tiles := []*image.RGBA{tile, tile, tile, tile}
	sheet := makeSheet(tiles)
	if sheet.Bounds().Dx() != 1290 || sheet.Bounds().Dy() != 1120 {
		t.Fatalf("sheet size %v", sheet.Bounds())
	}
	sheet = makeSheet(tiles[:3])
	if sheet.Bounds().Dy() != 560 {
		t.Fatalf("single-row sheet height %d", sheet.Bounds().Dy())
	}
	// Corner of the tile outside the thumbnail is grey; below the page is grey.
	if c := tile.RGBAAt(0, 559); c.R != 0xdd {
		t.Fatalf("tile background %v", c)
	}
}
