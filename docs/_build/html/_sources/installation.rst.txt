Installation
============

Tested on Linux Mint and Lubuntu (Debian-based syntax used here).

Create a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is highly recommended to use a *virtual environment* to install software with several requirements. We move to the folder where we will install the virtual environment (do not forget to replace the names between angle brackets). The created environment can be then activated using ``source``:

.. code-block:: python

   python3 -m venv <path_to_virtual_environment>/<venv_tiebenn>
   source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate


.. tip::

   You can add an alias by including the following line at the end of your ``~/.bashrc`` for quick access:

   .. code-block:: bash  

      alias <alias_name>='source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate'

   Save changes. Then reload:

   .. code-block:: bash

      exec bash

Installing TieBeNN and its Python dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you plan to work on a CPU machine, after activating the virtual environment by using the previously created alias, type:

.. code-block:: bash

   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu


If you are working on a GPU machine, then you can skip this line. Then, you can proceed to clone the TieBeNN repository and install it:

.. code-block:: bash

   git clone https://github.com/Cthuulhaa/tiebenn.git
   cd tiebenn

   pip install .
   pip install -r optional.txt

.. note::

   You are free to modify or comment-out those packages you will not need.

Installing NonLinLoc and setting paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*(Follow these steps in case NonLinLoc is not installed)* NonLinLoc must be individually compiled to make sure it is compatible with the machine where TieBeNN will be running. First we will create a convenient directory to install NonLinLoc and within it, we clone the NonLinLoc repository:

.. code-block:: bash

   mkdir -p <some_convenient_directory>
   cd <some_convenient_directory>

   git clone https://github.com/ut-beg-texnet/NonLinLoc.git

   cd NonLinLoc/src
   mkdir bin
   cmake .
   make
   cd bin
   cp Vel2Grid* Grid2* NLLoc ../../../

.. important::

   Do not use NonLinLoc's latest release directly, as it might contain unresolved bugs, whose fix are still unreleased.

Set NonLinLoc in your ``PATH``:

.. code-block:: bash

   echo 'export PATH=${PATH}:<some_convenient_directory>' >> ~/.bashrc
   exec bash

Install GMT
~~~~~~~~~~~

GMT version 6.0+ is recommended (officially 6.4+ for PyGMT). `The official GMT documentation <https://docs.generic-mapping-tools.org/dev/install.html>`_ has installation instructions, including instructions to migrate from earlier versions, and of course, a bunch of tutorials.
