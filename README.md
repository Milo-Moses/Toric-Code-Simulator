# Toric-Codes-Simulator

This project serves to give a simple introduction to the toric code as topological quantum computing, as well as to give a programatic environment in which to play with the code.

The bulk of the work on this project is mathematical. The literature on topological quantum computing can be very daunting, and sometimes it can be difficult to get what the big ideas are. The majority of time spent on this project was spent digesting all of that knowledge into a begginger-friendly format. Additionally, proofs to folklore and implitely assumed elementary results are formalized. The mathematical write up for the project can be found in the LaTeX Write Up folder.

The programatic implementation of the toric code is in Torus.py. While the final result is trivial, it was a huge task to try to find the correct was to implement the gates. Seeing as even small (in this case, 18) qubit systems take a large amount of memeory to store in classical computers, it took several tries to get the right format. We settle on using Qiskit. Despite not using the main QuantumCircuit feature of QisKit, we reference it here as simply a powerful linear algebra package.

The code is poorly commented. For documentation and an explination of how to use the code, see the JupyterNotebook, which walks step by step through all the key features.
