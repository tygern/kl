package main

import (
	"fmt"
	"html"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
)

type pdfInfo struct {
	pages  int
	title  *string
	author *string
}

// readInfo runs pdfinfo and parses Pages, Title and Author.
func readInfo(pdf string) (pdfInfo, error) {
	out, err := exec.Command("pdfinfo", "-enc", "UTF-8", pdf).Output()
	if err != nil {
		return pdfInfo{}, fmt.Errorf("pdfinfo failed on %s: %w", pdf, err)
	}
	return parseInfo(string(out))
}

func parseInfo(text string) (pdfInfo, error) {
	var info pdfInfo
	havePages := false
	for _, line := range strings.Split(text, "\n") {
		key, value, ok := strings.Cut(line, ":")
		if !ok {
			continue
		}
		value = strings.TrimSpace(value)
		switch strings.TrimSpace(key) {
		case "Pages":
			n, err := strconv.Atoi(value)
			if err != nil {
				return info, fmt.Errorf("bad pdfinfo Pages value %q", value)
			}
			info.pages, havePages = n, true
		case "Title":
			v := value
			info.title = &v
		case "Author":
			v := value
			info.author = &v
		}
	}
	if !havePages {
		return info, fmt.Errorf("pdfinfo output has no Pages line")
	}
	return info, nil
}

type wordBox struct {
	text                   string
	xMin, yMin, xMax, yMax float64
}

type pageBoxes struct {
	width, height float64
	words         []wordBox
}

var (
	pageTag = regexp.MustCompile(`<page width="([0-9.]+)" height="([0-9.]+)">`)
	wordTag = regexp.MustCompile(`<word xMin="([0-9.\-]+)" yMin="([0-9.\-]+)" xMax="([0-9.\-]+)" yMax="([0-9.\-]+)">([^<]*)</word>`)
)

// readWordBoxes runs pdftotext -bbox and returns the word boxes of each page.
func readWordBoxes(pdf string) ([]pageBoxes, error) {
	out, err := exec.Command("pdftotext", "-bbox", "-enc", "UTF-8", pdf, "-").Output()
	if err != nil {
		return nil, fmt.Errorf("pdftotext failed on %s: %w", pdf, err)
	}
	return parseBBox(string(out))
}

func parseBBox(text string) ([]pageBoxes, error) {
	var pages []pageBoxes
	chunks := strings.Split(text, "<page ")
	for _, chunk := range chunks[1:] {
		chunk = "<page " + chunk
		m := pageTag.FindStringSubmatch(chunk)
		if m == nil {
			return nil, fmt.Errorf("unparsable <page> tag in pdftotext output")
		}
		var p pageBoxes
		p.width, _ = strconv.ParseFloat(m[1], 64)
		p.height, _ = strconv.ParseFloat(m[2], 64)
		for _, w := range wordTag.FindAllStringSubmatch(chunk, -1) {
			var box wordBox
			var err error
			vals := [4]*float64{&box.xMin, &box.yMin, &box.xMax, &box.yMax}
			for i, v := range vals {
				if *v, err = strconv.ParseFloat(w[i+1], 64); err != nil {
					return nil, err
				}
			}
			box.text = html.UnescapeString(w[5])
			if strings.TrimSpace(box.text) == "" {
				continue
			}
			p.words = append(p.words, box)
		}
		pages = append(pages, p)
	}
	return pages, nil
}
