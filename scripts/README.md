# Instructions to launch a ParaView remote client/server connection

Execute launch script on Sherlock:  

`./home/groups/amarsden/karoline/launch_remote_paraview.sh`  

If you need more RAM, want to specify the duration of the job, modify the port number etc, you can modify the script accordingly.    
As soon as the job is running, your output will look like this:  

```
Waiting for client...
Connection URL: cs://<compute_node>:12345
Accepting connection(s): <compute_node>:12345
```

Copy the name of the compute node and open a new **local** terminal:  

`ssh -Y -L 12345:<compute_node>:12345 <username>@login.sherlock.stanford.edu`

Open your **local** ParaView installation (must be the exact same version as on Sherlock!) and connect to the server by finding the `Connect` button on the top ribbon:  

```
Add server
Name: Sherlock (or any you prefer)
Server type: Client / Server
Host: localhost
Port: 12345
Configure
Save
Connect
```

You should now see the successful connection in the Pipeline Browser and have access to your remote data. 
