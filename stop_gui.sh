#!/usr/bin/env bash
# Stop any running StRym GUI server by its listening port.
if command -v ss >/dev/null 2>&1; then
  pid=$(ss -ltnp 2>/dev/null | awk '/:3838 /{print $NF}' | grep -o 'pid=[0-9]*' | cut -d= -f2 | sort -u)
else
  pid=$(lsof -tiTCP:3838 -sTCP:LISTEN 2>/dev/null || true)
fi
[ -n "$pid" ] && kill $pid && echo "stopped $pid" || echo "no server running"
