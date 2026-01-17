package main

import (
	"errors"
	"fmt"
	"os"
	"strings"

	"scbamop/internal/app"
	"scbamop/internal/cli"
)

func main() {
	command, parseErr := cli.ParseArgs(os.Args[1:])
	if parseErr != nil {
		handleParseError(parseErr)
		return
	}

	switch command.Name {
	case "split":
		if command.Split == nil {
			fmt.Fprintln(os.Stderr, "Error: split command missing configuration")
			os.Exit(1)
		}
		runErr := app.RunSplit(*command.Split)
		if runErr != nil {
			fmt.Fprintln(os.Stderr, runErr)
			os.Exit(1)
		}
	default:
		cli.PrintGlobalUsage(os.Stderr)
		os.Exit(1)
	}
}

func handleParseError(parseErr error) {
	switch {
	case errors.Is(parseErr, cli.ErrShowUsage):
		cli.PrintGlobalUsage(os.Stderr)
		return
	case errors.Is(parseErr, cli.ErrShowSplitUsage):
		cli.PrintSplitUsage(os.Stderr)
		return
	}

	fmt.Fprintln(os.Stderr, "Error:", parseErr)
	if strings.HasPrefix(parseErr.Error(), "unknown command") {
		cli.PrintGlobalUsage(os.Stderr)
		os.Exit(1)
	}
	if strings.Contains(parseErr.Error(), "missing required arguments") {
		cli.PrintSplitUsage(os.Stderr)
	}
	os.Exit(1)
}
