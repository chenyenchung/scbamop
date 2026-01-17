package config

import (
	"strings"

	"scbamop/internal/logging"
)

type TagLocation int

const (
	TagLocationReadTag TagLocation = iota
	TagLocationReadName
)

type TagKind int

const (
	TagKindCellBarcode TagKind = iota
	TagKindUMI
)

type TagMeta struct {
	Location  TagLocation
	TagName   string
	Separator string
	Field     int
	Length    int
}

type SplitConfig struct {
	InputPath     string
	MetaPath      string
	OutputPrefix  string
	MapQThreshold int
	Deduplicate   bool
	DryRun        bool
	Verbose       bool
	LogLevel      logging.Level
	CellBarcode   TagMeta
	UMI           TagMeta
}

func DefaultTagMeta(tagName string) TagMeta {
	return TagMeta{
		Location:  TagLocationReadTag,
		TagName:   tagName,
		Separator: ",",
		Field:     1,
		Length:    20,
	}
}

func ApplyPlatform(tagMeta *TagMeta, platform string, kind TagKind) {
	if tagMeta == nil || platform == "" {
		return
	}
	platformLower := strings.ToLower(platform)
	switch kind {
	case TagKindCellBarcode:
		switch platformLower {
		case "10xv2":
			tagMeta.Length = 18
		case "scirnaseq3":
			tagMeta.Location = TagLocationReadName
			tagMeta.Length = 20
			tagMeta.Field = 1
		default:
			tagMeta.Length = 18
		}
	case TagKindUMI:
		switch platformLower {
		case "10xv2":
			tagMeta.TagName = "UB"
			tagMeta.Length = 10
		case "scirnaseq3":
			tagMeta.Location = TagLocationReadName
			tagMeta.Length = 8
			tagMeta.Field = 2
		default:
			tagMeta.TagName = "UB"
			tagMeta.Length = 12
		}
	}
}
