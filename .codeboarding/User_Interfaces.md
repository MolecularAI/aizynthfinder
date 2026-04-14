```mermaid

graph LR

    User_Interfaces["User Interfaces"]

    CLI_Interface["CLI Interface"]

    Web_GUI_Application["Web GUI Application"]

    GUI_Visualization_Modules["GUI Visualization Modules"]

    Core_AiZynthFinder["Core AiZynthFinder"]

    Analysis_Module["Analysis Module"]

    Context_Management["Context Management"]

    Utility_Module["Utility Module"]

    User_Interfaces -- "specializes" --> CLI_Interface

    User_Interfaces -- "specializes" --> Web_GUI_Application

    CLI_Interface -- "uses" --> Core_AiZynthFinder

    CLI_Interface -- "uses" --> Utility_Module

    Web_GUI_Application -- "uses" --> Core_AiZynthFinder

    Web_GUI_Application -- "orchestrates display of" --> GUI_Visualization_Modules

    Web_GUI_Application -- "retrieves results from" --> Analysis_Module

    Web_GUI_Application -- "uses" --> Utility_Module

    GUI_Visualization_Modules -- "visualizes data from" --> Analysis_Module

    GUI_Visualization_Modules -- "uses" --> Utility_Module

    Core_AiZynthFinder -- "depends on" --> Context_Management

    Core_AiZynthFinder -- "generates results for" --> Analysis_Module

    Core_AiZynthFinder -- "uses" --> Utility_Module

    Analysis_Module -- "uses" --> Utility_Module

    Context_Management -- "uses" --> Utility_Module

    click User_Interfaces href "https://github.com/MolecularAI/aizynthfinder/blob/master/.codeboarding//User_Interfaces.md" "Details"

    click Context_Management href "https://github.com/MolecularAI/aizynthfinder/blob/master/.codeboarding//Context_Management.md" "Details"

```



[![CodeBoarding](https://img.shields.io/badge/Generated%20by-CodeBoarding-9cf?style=flat-square)](https://github.com/CodeBoarding/GeneratedOnBoardings)[![Demo](https://img.shields.io/badge/Try%20our-Demo-blue?style=flat-square)](https://www.codeboarding.org/demo)[![Contact](https://img.shields.io/badge/Contact%20us%20-%20contact@codeboarding.org-lightgrey?style=flat-square)](mailto:contact@codeboarding.org)



## Details



The User Interfaces component in AiZynthFinder serves as the primary interaction layer, offering both a command-line interface (CLI) for scripting and a web-based graphical user interface (GUI) for interactive exploration. This design aligns with the Command Query Responsibility Segregation (CQRS) pattern, where the CLI handles "commands" (initiating searches) and the GUI focuses on "queries" (displaying and analyzing results).



### User Interfaces [[Expand]](./User_Interfaces.md)

The overarching component providing various means for users to interact with the AiZynthFinder application. It abstracts the different interaction modes.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.interfaces`</a>





### CLI Interface

Provides a command-line interface for batch processing, scripting, and direct execution of retrosynthetic searches. It parses arguments and orchestrates the core application logic.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/interfaces/aizynthcli.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.interfaces.aizynthcli`</a>





### Web GUI Application

Implements the web-based graphical user interface, managing user sessions, interactive search execution, and data presentation. It acts as the backend for the interactive user experience.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/interfaces/aizynthapp.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.interfaces.aizynthapp.AiZynthApp`</a>





### GUI Visualization Modules

Contains specific modules and utilities for rendering interactive visualizations within the web GUI, such as Pareto fronts, clustered results, and search tree displays.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.interfaces.gui`</a>





### Core AiZynthFinder

Encapsulates the primary retrosynthetic planning logic, coordinating search algorithms, policy application, and result generation. It's the central engine that both CLI and GUI interact with.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.aizynthfinder.AiZynthFinder`</a>





### Analysis Module

Responsible for post-processing and structuring the raw results from the search algorithms into meaningful data representations, such as route collections and tree analyses, suitable for display or further computation.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.analysis`</a>





### Context Management [[Expand]](./Context_Management.md)

Provides a centralized mechanism for managing application configuration, expansion policies, filter policies, scoring functions, and stock information. It acts as a dependency injection container for various strategies.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.context`</a>





### Utility Module

A collection of general-purpose helper functions and common utilities used across different parts of the application, including file handling, logging, and data manipulation.





**Related Classes/Methods**:



- <a href="https://github.com/MolecularAI/aizynthfinder/blob/master/aizynthfinder/aizynthfinder.py" target="_blank" rel="noopener noreferrer">`aizynthfinder.utils`</a>









### [FAQ](https://github.com/CodeBoarding/GeneratedOnBoardings/tree/main?tab=readme-ov-file#faq)