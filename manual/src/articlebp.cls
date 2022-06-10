\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{articlebp}[2016/05/23]

%% Load base class
\LoadClass{article}

%% Load packages
\usepackage[toc,page]{appendix}
\usepackage{geometry}
\usepackage{import}
\usepackage{titlesec}
\usepackage{cmap}
\usepackage{amsmath,amsfonts,amssymb}
\usepackage{varwidth}
\usepackage{graphicx}
\usepackage{setspace}
\usepackage[version=3]{mhchem}
\usepackage{graphicx}
\usepackage{latexsym}
\usepackage{keyval}
\usepackage{ifthen}
\usepackage{moreverb}
\usepackage{xcolor}
\usepackage{multicol}

%% Bib and Referencing
\usepackage[backend=bibtex,
bibstyle=numeric,
hyperref=true,
sorting=none]{biblatex}
\bibliography{bib/activity.bib}     %% CHANGE BIB FILE HERE
\usepackage{notoccite}

%% Page settings
\pagenumbering{arabic}
\geometry{a4paper,total={210mm,297mm},left=20mm,right=20mm,top=20mm,bottom=20mm}

%% Margin settings
\def\changemargin#1#2{\list{}{\rightmargin#2\leftmargin#1}\item[]}
\let\endchangemargin=\endlist

%% Paragraph settings
\usepackage{setspace}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.3ex plus 0.5ex minus 0.3ex}
\onehalfspacing
%%\doublespacing

%% Tikz package
\usepackage{tikz}
\usetikzlibrary{shapes,arrows}

%% Colours
\usepackage{color}
\definecolor{red}{rgb}{1.0,0,0}
\definecolor{green}{rgb}{0,1.0,0}
\definecolor{blue}{rgb}{0,0,1.0}
\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}

%% Listings for code
\usepackage{listings}
\lstloadlanguages{Fortran,Bash}
\lstset{frame=tb,
  aboveskip=3mm,
  belowskip=3mm,
  showstringspaces=false,
  columns=flexible,
  basicstyle={\small\ttfamily}
}
%% Listing styles
\lstdefinestyle{sBash}
{
  numbers=left,
  stepnumber=1,
  language=bash,
  tabsize=4,
  keywordstyle=\color{red}
}
\lstdefinestyle{sFortran}
{
  numbers=left,
  stepnumber=1,
  language=Fortran,
  tabsize=4,
  numberstyle=\tiny\color{gray},
  keywordstyle=\color{blue},
  commentstyle=\color{dkgreen},
  stringstyle=\color{mauve}
}

%% Algorithm + Pseduocode
\usepackage{algorithmicx}
\usepackage{algpseudocode}
\algdef{SE}[DOWHILE]{Do}{doWhile}{\algorithmicdo}[1]{\algorithmicwhile\ #1}%


%%------------------------------------------------------------------------------
%%------------------------------------------------------------------------------
%% Define new commands

%% Thesis
%%-----------------------------
\def\changemargin#1#2{\list{}{\rightmargin#2\leftmargin#1}\item[]}
\let\endchangemargin=\endlist
\newcommand{\abstractpreamble}[1]{\textit{#1}}
\newcommand{\subheadings}{\noindent\textbf}


%% Articles
%%-----------------------------
\newcommand{\R}{\mathbb{R}}



%% CVs
%%-----------------------------
\newcommand{\cvsection}[1]{\large{\bf{#1}}}
\newcommand{\cvsubsection}[1]{\normalsize{\bf{#1}}}
\newcommand{\cvtext}[1]{\normalsize{#1}}











%%------------------------------------------------------------------------------
%%------------------------------------------------------------------------------
