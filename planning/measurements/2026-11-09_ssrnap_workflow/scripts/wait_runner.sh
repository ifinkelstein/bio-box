#!/bin/sh
# Local helper: wait until the remote runner log has its final "exit" line, then print it.
until ssh -o BatchMode=yes if2453@wilkcomp02.ccbb.utexas.edu "grep -q '^exit' /tmp/if2453/ssrnap/logs/runner.log"; do
  sleep 20
done
ssh -o BatchMode=yes if2453@wilkcomp02.ccbb.utexas.edu "cat /tmp/if2453/ssrnap/logs/runner.log"
