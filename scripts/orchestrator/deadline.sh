#!/bin/sh
# Wake the session after N minutes with a label, so a silent job or agent gets checked.
sleep $(( $1 * 60 ))
echo "DEADLINE reached after $1 min: $2"
