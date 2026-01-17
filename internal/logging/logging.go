package logging

import (
	"fmt"
	"os"
)

type Level int

const (
	Error   Level = 0
	Warning Level = 1
	Info    Level = 3
	Debug   Level = 5
)

var levelNames = []string{"ERROR", "WARNING", "", "INFO", "", "DEBUG"}

type Logger struct {
	Level Level
}

func (logger Logger) Logf(level Level, format string, args ...interface{}) {
	if level > logger.Level {
		return
	}
	levelName := ""
	if int(level) >= 0 && int(level) < len(levelNames) {
		levelName = levelNames[level]
	}
	fmt.Fprintf(os.Stderr, "[%s] ", levelName)
	fmt.Fprintf(os.Stderr, format, args...)
	fmt.Fprintln(os.Stderr)
}
