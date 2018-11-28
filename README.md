# Quantum Control with a Genetic Algorithm: An Example

This repository is intended as an example for researchers learning how to use a genetic algorithm (GA). I'll describe the basics of a GA, our specificic problem, and how we applied the GA to our problem.

Our code was first developed in C++ in simulation and then translated to C# for the actual experimental implementation. That code is available here but since it might not be useful for everyone, we've also provided a pseudo-code version.

# Controlling an Electron's "Path" to Ionization

I'll try to summarize the problem without too many technical details. If you're interested in those details, you can find a full report of our experiments published in the journal Physical Review A:

* Vincent C. Gregoric, Xinyue Kang, Zhimin Cheryl Liu, Zoe A. Rowley, Thomas J.
Carroll, and Michael W. Noel. Quantum control via a genetic algorithm of the field
ionization pathway of a Rydberg electron. Phys. Rev. A, 96(2):023403, August 2017.
* Vincent C. Gregoric, Jason J. Bennett, Bianca R. Gualtieri, Ankitha Kan-
nad, Zhimin Cheryl Liu, Zoe A. Rowley, Thomas J. Carroll, and Michael W.
Noel. Improving the state selectivity of field ionization with quantum control.
arXiv:1806.00889 [physics], June 2018. Accepted to Phys. Rev. A.

Imagine that you would like to measure the energy of the outer-most electron in a "highly excited" atom called a Rydberg atom. Don't forget that in quantum mechanics the possible results of this measurement can only have particular values! What if two of the possible energies that the electron could have are very close together? It might be hard to distinguish between them with your measurement.

This is, in fact, what can happen and it's the problem we aim to solve with our GA. A fairly standard technique for measuring the energy of the electron is to rip it off of the atom with a big pulse of electric field. Imagine that the electron is stuck in the "gravity well" of the atom. The farther down the well, the more energy it will take to rip it away from the atom. It takes us a longer time to give it that energy, so electrons farther down the well will take longer to arrive at our detector. In this way, we convert a difference in energy (deeper or shallower in the well) to a difference in arrival time at our detector (later or earlier).

Of course, it's quite a bit more complicated than that simple description. In reality, the electron follows a "path" through many hundreds or thousands of energy levels on it's way out of the well. This being quantum mechanics, there are many opportunities for "interference."  To make a long story short (too late!), if we could take that big pulse of electric field and put some wiggles on top of it we should be able to control what path the electron takes.

So our solution is to use a GA to evolve the wiggles so that electrons that start in close-together energy levels take diverging paths and end up arriving at our detector at different times -- more different than they would without our evolved wiggles.

# Why a Genetic Algorithm?

There are lots of optimization methods out there. We had some good reasons for choosing a GA. But it's also true that we tried it and it worked, so we stuck with it! But the reasons are important and might help you to decide whether a GA is what you want for your problem.

A GA is a "stochastic" optimization method. That basically means that there is some randomness used to generate possible solutions. This is great if you don't know what your solution should look like. The alogrithm does a pretty good job searching without your guidance.

The GA also automatically takes into account lots of variables that you might not be able to characterize very well. In our case, for example, if we wanted to calculate our solution from first principles we would need to know very well the value of our inhomogeneous magnetic field. We don't know it and it's really hard to measure. So we just let the GA evolve around it!

