#!/bin/zsh
# watch_machine.sh -- emit an event when the other agent's heavy stage finishes,
# and when swap enters the range that killed both Python processes last time.
#
# Two agents share this 36 GB box. The ALT build must not start tetgen or mmg
# while PREFERRED's mmg is running. Rather than poll by hand, emit one line when
# it exits (that is the go signal) and a line whenever swap is critical.
# Coverage matters more than tidiness here: a filter that only reported success
# would stay silent through a crash, which looks identical to "still running".
while true; do
  sw=$(sysctl -n vm.swapusage | sed 's/.*used = \([0-9.]*\)M.*/\1/')
  st=$(sysctl -n vm.swapusage | sed 's/.*total = \([0-9.]*\)M.*/\1/')
  pct=$(echo "$sw $st" | awk '{printf "%d", 100*$1/$2}')
  mm=$(ps -Ao rss=,command= | grep mmg3d | grep -v grep | awk '{printf "%.1f", $1/1048576}' | head -1)
  tg=$(ps -Ao rss=,command= | grep -E "tetgen|fill_domain" | grep -v grep | awk '{printf "%.1f", $1/1048576}' | head -1)
  if [ -z "$mm" ] && [ -z "$tg" ]; then
    echo "HEAVY-CLEAR: no mmg/tetgen running; swap ${pct}% (${sw}M/${st}M) -- ALT heavy stages may start"
    exit 0
  fi
  if [ "$pct" -ge 96 ]; then
    echo "SWAP-CRITICAL ${pct}% (${sw}M/${st}M) mmg=${mm:-none}GB tetgen=${tg:-none}GB -- crash risk, do not add load"
  fi
  sleep 60
done
