# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Sunday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
10/02 and 11/02: 
Branch 'testing' created because no permission to edit master branch.
1) Lennard Johnes function for calculating the potential between a pair of molecules created. Two versions created, one with handling errors and exceptions and one without (version for reading fast what's going on in the code)
2) file main_module created, contains constructors of class Molecule (describing a single molecule/atom) and Gas describing collection of instances of Molecule class. Methods added to the classes, eg to obtain position, velocity of molecules and forces experienced due to Lennard Johnes potential
3) simulation.py file created. This is the simulation where we create molecules and enviroment (eg can specify temperature, pressure etc). Other files are imported and used by this one and run 'behind the scenes'.
11/02-13/02:
1) code expanded to support 2D
2) time evolution of the system added
3) periodic boundary conditions added
4) calculation of total energy, potential energy and kinetic energy and check of conservation of energy added

(due 14 February 2021, 23:59)


## Week 2
(due 21 February 2021, 23:59)


## Week 3
(due 28 February 2021, 23:59)


## Week 4
(due 7 March 2021, 23:59)


## Week 5
(due 14 March 2021, 23:59)
