**CONTAINER USAGE**

With [Singularity](https://sylabs.io/docs/), you will want to pull the latest container version and then run the container:

    singularity pull library://sestaton/default/tephra
    singularity run --bind $PWD/db:/db tephra_latest.sif

That assumes you have your data files for analysis (genome, config file, etc.) in a directory called `db` in the current working directory. After that, change to the `/db` directory and run your analysis.

---

With [Docker](https://www.docker.com/), you can create a container to run Tephra with the following command:

    docker run -it --name tephra-con -v $(pwd)/db:/db:Z sestaton/tephra

That will create a container called `tephra-con` and start an interactive shell. The above assumes you have a directory called `db` in the working directory that contains your database files and the Tephra configuration. To run the full analysis, change to the mounted directory with `cd /db` in your container and run the following command:

    tephra all -c tephra_config.yml

I recommend using `nohup` and then logging out, which will allow you to leave the container running in the background.

If you cannot use Singularity or Docker, please see the [INSTALL](https://github.com/sestaton/tephra/blob/master/INSTALL.md) file included with this distribution to install Tephra on various operating systems.

**SUPPORT AND DOCUMENTATION**

You can get usage information at the command line with the following command:

    perldoc tephra

The `tephra` program will also print a diagnostic help message when executed with no arguments, and display the available subcommands.

You can also look for information at:

    Tephra wiki
            https://github.com/sestaton/tephra/wiki

    Tephra issue tracker
            https://github.com/sestaton/tephra/issues

 

