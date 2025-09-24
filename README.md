# A WORK IN PROGRESS

Its been a while since I first put this together and realized there were some clear problems. I intend to revamp this resource in the coming months, but for now, I put together a much simplified example for git here:

[github.com/nsantantonio/exampleGit](https://github.com/nsantantonio/exampleGit)


# Reproducible Science

Welcome to my take on how one might use git and github as a way to (re)produce science. Here I will briefly present how to use several tools to build reproducible science projects. 

# Gold Standard

The gold standard for science, or a scientific manuscript rather, is for complete reproducability. I.e. another researcher should be able to reproduce your experiment exactly, given access to the correct starting materials. We know the outcomes will vary (vartiation due to sampling, environmental conditions, etc.), but the experimental design should be able to be replicated exactly. 

# A note on spreadsheets

Spreadsheets have their place. At one point in my career, I refused to use them for even the simplest tasks (I also used to use \LaTeX to make presentations, total waste of time). I now use them regularly, although I generally abstain from doing any data manipulation in them if I need that data for anything important. I encourage you to 

# Multiple collaborators

Git has awesome features for multiple people working on the same project. However, outside of multi-person software development (arguably what git was designed for), I have found little benefit to sharing git repositories outside of simply making everything publically available and having an official public version while I tinker locally.  

Internally, git can save you a lot of heartache through version control. Ever break a script and cant find the problem? You can either revert to an earlier version that worked or continue to tinker without fear that you will make things worse. 

For collaborative writing, I find Overleaf to be more useful than a git repo. Overleaf is a more "google docs" type platform for multiple editors to \LaTeX documents

## Basic git 101

	git init

	git add <filename>
		- alternatively -A for everything unadded (including untracked files), and -u for adding only for currently tracked files

	git commit -m "my message to know what changes were made"

	git push origin master

## Branches

I am going to avoid from discussion of branches here, but the idea is you can branch off, make a bunch of changes and then when you are satisfied with the result, you merge the branch back into the master. This also creates for a safe place for multiple people to be working on the same project at the same time. Then once everyone's various changes are shown to work together, then you merge the branch back to the master. I rarely do this unless I am working on large complicated projects (probably poor practice). See a list of tutorials below

## Some useful git commands

check what changes have been made but not added
	
	git status

If you want to know what has been changed since the last time you worked on the project,
	
	git diff <filename>

If you accidentally staged a file for a commit (i.e. added it)

	git reset <filename>

# For the Workshop

I have put a small example together to show how a git repository can be used to make reproducible science in the form of data analysis integrated into a manuscript. This is not th eonly use for git repositories, but I think it will be the most useful for this audience, i.e. Grad students in plant sciences.

The main files are 
* 

## Downloads and installs
Please make sure to download and/or install the following:

- git [https://git-scm.com/downloads]; or if you are on a linux machine, install from your package manager, e.g. apt (debian, ubuntu, etc) or yum (centos, etc.). If you are on a mac, homebrew I think is still most common for package management. I assume windows works similarly to unix machines (via windows powershell?) but I haven't tried it yet**. I intend to try before seminar tomorrow. 

- R [https://cran.r-project.org/]; this is just if you want to run any of the scripts
	- the devtools package is also very useful for interacting with github packages in R (as well as making your own packages!) see below in optional


Please also make an account on Github here:
- Github <https://github.com/>

# Optional:

- Rstudio <https://posit.co/download/rstudio-desktop/>; I will not be using Rstudio, but I find many are familiar with it, so you may want to have it on your machine

- devtools is an R library that is useful for making packages and installing packages from github. Install with install.packages("devtools"). Note, you may have issues with install as devtools is a huge package with a lot of dependencies. I.e. install at your own risk!

- \LaTeX; If you really want to get into the weeds, also install \LaTeX. \LaTeX is a typesetting language used for making very nice documents, (most of your textbooks were probably made using some flavor of \TeX). \LaTeX is awesome for reproduceable science because you can I recommend installing the latest version of TeX Live <https://www.tug.org/texlive/>, but you can often also install from a package manager. 


*Note, if you cannot get all of this to work on your system, you can still attend and listen in. I do believe having a machine in front of you and trying to follow along may be useful for some.

**it appears there are a lot of GUIs (graphic user interfaces) for git on windows, and these may facilitate your use of git, but they seem unnecessary to me beyond familiarity of using a mouse.

# Some links to tutorials

I haven't vetted these super well, but read throught them briefly. They give you the basics. 

Remember, you can google (or LLM) just about anything these days. Just ask your favorite search engine / AI, it can probably lead you to the right solution. 

1. Official git documentation <https://git-scm.com/docs/gittutorial>, surprisingly useful, although may be tough for very early users. 
2. <https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners>
3. <http://www.w3schools.com/GIT/default.asp?remote=github>

# Other useful links 

- The summaryTools package [github.com/nsantantonio/summaryTools] has a couple functions from making latex tables in R, names `xtable2()` and `longtable()`. You can write dataframes in R to latex tables, so that if you change your analysis, you can automatically update tables in your manuscript.
