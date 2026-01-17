package cli

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
	"unicode"

	"scbamop/internal/config"
	"scbamop/internal/logging"
)

var ErrShowUsage = errors.New("show usage")
var ErrShowSplitUsage = errors.New("show split usage")

type Command struct {
	Name  string
	Split *config.SplitConfig
}

func ParseArgs(args []string) (Command, error) {
	if len(args) == 0 {
		return Command{}, ErrShowUsage
	}

	subcommand := args[0]
	if subcommand == "--help" || subcommand == "-h" {
		return Command{}, ErrShowUsage
	}

	if subcommand == "split" {
		splitConfig, parseErr := ParseSplitArgs(args[1:])
		if parseErr != nil {
			return Command{}, parseErr
		}
		return Command{Name: "split", Split: &splitConfig}, nil
	}

	return Command{}, fmt.Errorf("unknown command %q", subcommand)
}

func ParseSplitArgs(args []string) (config.SplitConfig, error) {
	splitConfig := config.SplitConfig{
		OutputPrefix: "./",
		LogLevel:     logging.Warning,
		CellBarcode:  config.DefaultTagMeta("CB"),
		UMI:          config.DefaultTagMeta("UB"),
	}

	argumentIndex := 0
	for argumentIndex < len(args) {
		argument := args[argumentIndex]
		switch {
		case argument == "-f" || argument == "--file":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.InputPath = value
		case argument == "-m" || argument == "--meta":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.MetaPath = value
		case argument == "-o" || argument == "--output":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.OutputPrefix = value
		case argument == "-q" || argument == "--mapq":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			mapqValue, parseErr := parseNonNegativeInt(value, "MAPQ threshold")
			if parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.MapQThreshold = mapqValue
		case argument == "-p" || argument == "--platform":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			config.ApplyPlatform(&splitConfig.CellBarcode, value, config.TagKindCellBarcode)
			config.ApplyPlatform(&splitConfig.UMI, value, config.TagKindUMI)
		case argument == "-d" || argument == "--dedup":
			splitConfig.Deduplicate = true
		case argument == "-b" || argument == "--cbc-location":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			parseErr = applyTagLocation(value, &splitConfig.CellBarcode)
			if parseErr != nil {
				return splitConfig, parseErr
			}
		case argument == "-L" || argument == "--cbc-length":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			lengthValue, parseErr := parsePositiveLength(value, "cell barcode")
			if parseErr != nil {
				return splitConfig, parseErr
			}
			if parseErr := confirmLargeLength(lengthValue, "Cell barcode"); parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.CellBarcode.Length = lengthValue
		case argument == "-u" || argument == "--umi-location":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			parseErr = applyTagLocation(value, &splitConfig.UMI)
			if parseErr != nil {
				return splitConfig, parseErr
			}
		case argument == "-l" || argument == "--umi-length":
			value, parseErr := nextArgument(args, &argumentIndex, argument)
			if parseErr != nil {
				return splitConfig, parseErr
			}
			lengthValue, parseErr := parsePositiveLength(value, "UMI")
			if parseErr != nil {
				return splitConfig, parseErr
			}
			if parseErr := confirmLargeLength(lengthValue, "UMI"); parseErr != nil {
				return splitConfig, parseErr
			}
			splitConfig.UMI.Length = lengthValue
		case argument == "-n" || argument == "--dry-run":
			splitConfig.DryRun = true
		case argument == "-v" || argument == "--verbose":
			levelValue := ""
			if argumentIndex+1 < len(args) && isSingleDigit(args[argumentIndex+1]) {
				levelValue = args[argumentIndex+1]
				argumentIndex++
			}
			parseErr := applyVerbose(&splitConfig, levelValue)
			if parseErr != nil {
				return splitConfig, parseErr
			}
		case argument == "-h" || argument == "--help":
			return splitConfig, ErrShowSplitUsage
		case strings.HasPrefix(argument, "--verbose="):
			levelValue := strings.TrimPrefix(argument, "--verbose=")
			parseErr := applyVerbose(&splitConfig, levelValue)
			if parseErr != nil {
				return splitConfig, parseErr
			}
		case strings.HasPrefix(argument, "-v") && len(argument) > 2:
			levelValue := strings.TrimPrefix(argument, "-v")
			parseErr := applyVerbose(&splitConfig, levelValue)
			if parseErr != nil {
				return splitConfig, parseErr
			}
		default:
			return splitConfig, fmt.Errorf("unknown option %q", argument)
		}
		argumentIndex++
	}

	if splitConfig.InputPath == "" || splitConfig.MetaPath == "" {
		return splitConfig, fmt.Errorf("missing required arguments (-f and -m)")
	}

	if splitConfig.OutputPrefix == "" {
		splitConfig.OutputPrefix = "./"
	}
	if !strings.HasSuffix(splitConfig.OutputPrefix, "/") {
		splitConfig.OutputPrefix = splitConfig.OutputPrefix + "/"
	}

	if splitConfig.CellBarcode.Length <= 0 {
		return splitConfig, fmt.Errorf("cell barcode length must be larger than 0")
	}
	if splitConfig.UMI.Length <= 0 {
		return splitConfig, fmt.Errorf("UMI length must be larger than 0")
	}

	return splitConfig, nil
}

func nextArgument(args []string, argumentIndex *int, option string) (string, error) {
	nextIndex := *argumentIndex + 1
	if nextIndex >= len(args) {
		return "", fmt.Errorf("option %s requires an argument", option)
	}
	*argumentIndex = nextIndex
	return args[nextIndex], nil
}

func parseNonNegativeInt(value string, label string) (int, error) {
	parsedValue, parseErr := strconv.ParseInt(value, 10, 64)
	if parseErr != nil || parsedValue < 0 {
		return 0, fmt.Errorf("invalid %s: %s", label, value)
	}
	return int(parsedValue), nil
}

func parsePositiveLength(value string, label string) (int, error) {
	parsedValue, parseErr := strconv.ParseInt(value, 10, 64)
	if parseErr != nil || parsedValue <= 0 {
		return 0, fmt.Errorf("%s length must be larger than 0", label)
	}
	return int(parsedValue), nil
}

func applyTagLocation(value string, tagMeta *config.TagMeta) error {
	parsedValue, parseErr := strconv.ParseInt(value, 10, 64)
	if parseErr != nil {
		return setTagName(value, tagMeta)
	}
	if parsedValue < 0 {
		return fmt.Errorf("invalid tag field: %s", value)
	}
	if parsedValue == 0 {
		return setTagName(value, tagMeta)
	}
	tagMeta.Location = config.TagLocationReadName
	tagMeta.Field = int(parsedValue)
	return nil
}

func setTagName(value string, tagMeta *config.TagMeta) error {
	if len(value) > 2 {
		return fmt.Errorf("tag name too long (max 2 chars): %s", value)
	}
	tagMeta.TagName = value
	return nil
}

func applyVerbose(splitConfig *config.SplitConfig, levelValue string) error {
	splitConfig.Verbose = true
	if levelValue == "" {
		splitConfig.LogLevel = logging.Info
		return nil
	}
	parsedValue, parseErr := strconv.ParseInt(levelValue, 10, 64)
	if parseErr != nil {
		return fmt.Errorf("invalid verbosity level: %s", levelValue)
	}
	switch {
	case parsedValue == 0:
		splitConfig.LogLevel = logging.Info
	case parsedValue > 0 && parsedValue < 5:
		splitConfig.LogLevel = logging.Level(parsedValue)
	default:
		splitConfig.LogLevel = logging.Debug
	}
	return nil
}

func isSingleDigit(value string) bool {
	if len(value) != 1 {
		return false
	}
	return unicode.IsDigit([]rune(value)[0])
}

func confirmLargeLength(length int, label string) error {
	const maxRecommendedLength = 31
	if length <= maxRecommendedLength {
		return nil
	}
	if !isTerminal(os.Stdin) {
		fmt.Fprintf(os.Stderr, "Warning: %s length %d exceeds %d. Proceeding without prompt.\n", label, length, maxRecommendedLength)
		return nil
	}

	reader := bufio.NewReader(os.Stdin)
	for {
		fmt.Fprintf(os.Stderr, "Warning: %s length %d exceeds %d. Continue? [y/n] ", label, length, maxRecommendedLength)
		input, err := reader.ReadString('\n')
		if err != nil {
			return fmt.Errorf("failed to read confirmation: %w", err)
		}
		answer := strings.TrimSpace(strings.ToLower(input))
		switch answer {
		case "y", "yes":
			return nil
		case "n", "no":
			return fmt.Errorf("aborted due to %s length %d", strings.ToLower(label), length)
		default:
			fmt.Fprintln(os.Stderr, "Please answer 'y' or 'n'.")
		}
	}
}

func isTerminal(file *os.File) bool {
	info, err := file.Stat()
	if err != nil {
		return false
	}
	return info.Mode()&os.ModeCharDevice != 0
}
