#!/usr/bin/env bash
# Stop the StRym GUI server listening on STRYM_PORT (default 3838).
set -euo pipefail

STRYM_PORT="${STRYM_PORT:-3838}"
if [[ ! "$STRYM_PORT" =~ ^[0-9]+$ ]] || (( ${#STRYM_PORT} > 5 )); then
  echo "Invalid STRYM_PORT: $STRYM_PORT. Use a number from 1 to 65535." >&2
  exit 2
fi
STRYM_PORT=$((10#$STRYM_PORT))
if (( STRYM_PORT < 1 || STRYM_PORT > 65535 )); then
  echo "Invalid STRYM_PORT: $STRYM_PORT. Use a number from 1 to 65535." >&2
  exit 2
fi

if command -v ss >/dev/null 2>&1; then
  pid=$(ss -ltnp 2>/dev/null | awk -v port=":${STRYM_PORT}" 'index($4, port) == length($4) - length(port) + 1 { print $NF }' | grep -o 'pid=[0-9]*' | cut -d= -f2 | sort -u || true)
else
  pid=$(lsof -tiTCP:"${STRYM_PORT}" -sTCP:LISTEN 2>/dev/null || true)
fi
[ -n "$pid" ] && kill $pid && echo "stopped $pid on port ${STRYM_PORT}" || echo "no server running on port ${STRYM_PORT}"
