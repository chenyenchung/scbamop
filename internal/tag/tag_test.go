package tag

import (
	"testing"

	"github.com/biogo/hts/sam"

	"scbamop/internal/config"
)

func TestExtractFromTag(t *testing.T) {
	aux, err := sam.NewAux(sam.NewTag("CB"), "CELLBARCODE")
	if err != nil {
		t.Fatalf("failed to build aux tag: %v", err)
	}
	record := sam.Record{Name: "read1", AuxFields: []sam.Aux{aux}}
	value, err := Extract(&record, config.TagMeta{Location: config.TagLocationReadTag, TagName: "CB", Length: 20})
	if err != nil {
		t.Fatalf("expected tag value, got %v", err)
	}
	if value != "CELLBARCODE" {
		t.Fatalf("unexpected tag value: %s", value)
	}
}

func TestValidateLength(t *testing.T) {
	if err := ValidateLength("ABCDE", 5, "cell barcode"); err != nil {
		t.Fatalf("expected length validation to succeed, got %v", err)
	}
	if err := ValidateLength("AB", 5, "cell barcode"); err == nil {
		t.Fatalf("expected length validation error")
	}
}

func TestExtractFromReadName(t *testing.T) {
	record := sam.Record{Name: "foo|BARCODE|baz"}
	value, err := Extract(&record, config.TagMeta{Location: config.TagLocationReadName, Separator: "|", Field: 2, Length: 20})
	if err != nil {
		t.Fatalf("expected read name value, got %v", err)
	}
	if value != "BARCODE" {
		t.Fatalf("unexpected read name value: %s", value)
	}
}

func TestExtractMissingTag(t *testing.T) {
	record := sam.Record{Name: "read3"}
	_, err := Extract(&record, config.TagMeta{Location: config.TagLocationReadTag, TagName: "CB", Length: 20})
	if err == nil {
		t.Fatalf("expected missing tag error")
	}
}

func TestExtractUnsupportedTagType(t *testing.T) {
	aux, err := sam.NewAux(sam.NewTag("CB"), int32(42))
	if err != nil {
		t.Fatalf("failed to build aux tag: %v", err)
	}
	record := sam.Record{Name: "read4", AuxFields: []sam.Aux{aux}}
	_, err = Extract(&record, config.TagMeta{Location: config.TagLocationReadTag, TagName: "CB", Length: 20})
	if err == nil {
		t.Fatalf("expected unsupported tag type error")
	}
}

func TestExtractReadNameFieldMissing(t *testing.T) {
	record := sam.Record{Name: "foo|bar"}
	_, err := Extract(&record, config.TagMeta{Location: config.TagLocationReadName, Separator: "|", Field: 3, Length: 20})
	if err == nil {
		t.Fatalf("expected missing field error")
	}
}
