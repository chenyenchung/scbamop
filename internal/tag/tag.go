package tag

import (
	"fmt"
	"strings"

	"github.com/biogo/hts/sam"

	"scbamop/internal/config"
)

func Extract(record *sam.Record, tagMeta config.TagMeta) (string, error) {
	if record == nil {
		return "", fmt.Errorf("missing record")
	}
	switch tagMeta.Location {
	case config.TagLocationReadTag:
		return extractFromTag(record, tagMeta)
	case config.TagLocationReadName:
		return extractFromReadName(record.Name, tagMeta)
	default:
		return "", fmt.Errorf("unknown tag location")
	}
}

func extractFromTag(record *sam.Record, tagMeta config.TagMeta) (string, error) {
	if len(tagMeta.TagName) != 2 {
		return "", fmt.Errorf("tag name must be 2 characters")
	}
	auxValue, ok := record.Tag([]byte(tagMeta.TagName))
	if !ok {
		return "", fmt.Errorf("tag not found")
	}
	value := auxValue.Value()
	var tagValue string
	switch typedValue := value.(type) {
	case string:
		tagValue = typedValue
	case []byte:
		tagValue = string(typedValue)
	default:
		return "", fmt.Errorf("unsupported tag type")
	}
	return tagValue, nil
}

func extractFromReadName(name string, tagMeta config.TagMeta) (string, error) {
	if name == "" {
		return "", fmt.Errorf("empty read name")
	}
	if tagMeta.Separator == "" {
		return "", fmt.Errorf("missing separator")
	}
	if name[0] == tagMeta.Separator[0] {
		return "", fmt.Errorf("read name starts with separator")
	}

	fields := strings.Split(name, tagMeta.Separator)
	if tagMeta.Field <= 0 || tagMeta.Field > len(fields) {
		return "", fmt.Errorf("field not found")
	}
	value := fields[tagMeta.Field-1]
	return value, nil
}

func ValidateLength(value string, expected int, label string) error {
	if expected <= 0 {
		return nil
	}
	if len(value) != expected {
		return fmt.Errorf("%s length %d does not match expected %d", label, len(value), expected)
	}
	return nil
}
