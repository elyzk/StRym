#!/usr/bin/env bash
# Stop any running TuxNet-R GUI server (by listening port, so it can never match itself).
pid=$(ss -ltnp 2>/dev/null | awk '/:3838 /{print $NF}' | grep -o 'pid=[0-9]*' | cut -d= -f2 | sort -u)
[ -n "$pid" ] && kill $pid && echo "stopped $pid" || echo "no server running"
