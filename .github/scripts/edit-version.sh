#!/bin/bash
set -euo pipefail

# This script updates the version in __init__.py only
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <version>"
  exit 1
fi

echo "Updating version to: $1"
sed -i "s/^__version__ = .*/__version__ = \"$1\"/" ./src/ipp/__init__.py
echo "Version updated successfully to $1 in __init__.py."
