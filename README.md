This is a self-contained client distribution of the NCATS scaffold
hopper. A brief description of the tool is available
[here](https://tripod.nih.gov/?p=483). To run the scaffold hopper
client locally in your environment, you'll need to have (i) [Java
1.7](https://java.com/download/) or higher installed and (ii)
associated scaffold hopper web services. If you don't want to setup
your own services locally, then you can use our public instance at
https://tripod.nih.gov/chembl-v23. Simply launch the scaffold hopper
client without arguments:

```
./run.sh
```

Note there are two arguments here.

If you prefer to have the services running locally (so that
structure searching doesn't go outside of your organization), then you
can use our [docker
image](https://hub.docker.com/r/ncats/hopper-server/). First, you need
to install a recent version of
[docker](https://www.docker.com/community-edition). Then run the
following command to pull down the image:

```
docker pull ncats/hopper-server
```

Now to start the server, simply issue the command:

```
docker run -p 8080:8080 -it ncats/hopper-server:chembl24 /hopper-server
```

This assumes that port 8080 is available in the host environment. The
server will require 10-15 minutes to initialize. Now to make use of
your local web services, simply point the client locally as follows:

```
./run.sh http://localhost:8080 chembl-v24
```

Please report problems and/or suggestions to nguyenda@mail.nih.gov.
