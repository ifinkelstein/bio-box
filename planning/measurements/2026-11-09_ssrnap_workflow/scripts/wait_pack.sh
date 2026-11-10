#!/bin/sh
# Local helper: wait until the remote pack log has its final "exit" line, then print the tail.
until ssh -o BatchMode=yes if2453@wilkcomp02.ccbb.utexas.edu "grep -q '^exit' /tmp/if2453/ssrnap/logs/pack2.log"; do
  sleep 20
done
ssh -o BatchMode=yes if2453@wilkcomp02.ccbb.utexas.edu "tail -40 /tmp/if2453/ssrnap/logs/pack2.log"
