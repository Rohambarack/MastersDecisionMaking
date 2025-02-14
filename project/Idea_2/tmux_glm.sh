#!/usr/bin/env bash

#install packages
source setup_brms.sh

SESSION_NAME="T"

# Check if the session already exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching to it."
    tmux attach-session -t $SESSION_NAME
else
    # Create a new session and panes 
    # (there must be a more automated way with a loop)
    tmux new-session -d -s $SESSION_NAME
    
    tmux split-window -h
    tmux split-window -h
    tmux select-pane -t 0
    tmux split-window -h
    tmux select-pane -t 0
    tmux split-window -v
    tmux select-pane -t 2
    tmux split-window -v
    tmux select-pane -t 4
    tmux split-window -v
    tmux select-pane -t 6
    tmux split-window -v

    
    #loop for seeds and send them to the script
    for i in {0..7}
    do
      #define seed
      SEED=$((i + 1))
      # Send a command to the pane, make R read SEED
      tmux send-keys -t $i 'Rscript hier_glm.R ' $SEED C-m
    
    done
 # Attach to the created session
  tmux attach-session -t $SESSION_NAME
fi

