#!/usr/bin/env bash
set -euo pipefail

mkdir -p java-bin
javac -d java-bin $(find java-src -name "*.java")
