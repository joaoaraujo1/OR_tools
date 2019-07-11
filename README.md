# Operational Research Tools
This repository features a library of functions I have created when doing my Minor in Statistics and OR and implements some typical algorithms and methods used in operational research and multivariate exploratory analyses. All code is written in MATLAB

## Optimization - Flows and Networks (IP)
- Pape, D'Esopo and Moore algorithm (PDM): Shortest or longest path between a node and every other node on a graph (also detects cycles)
- [Floyd-Warshall](https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm) algorithm: Shortest or longest path any pair of nodes on a graph (also detects cycles)
- [Ford-Fulkerson](https://en.wikipedia.org/wiki/Ford%E2%80%93Fulkerson_algorithm): Maximum feasable and conservative flow in a network
- [Out-of-Kilter](https://en.wikipedia.org/wiki/Out-of-kilter_algorithm): Minimum cost flow in a network
- Greedy1: Lower limit value for the [Knapsack Problem](https://en.wikipedia.org/wiki/Knapsack_problem)
- [Facility Location Problem](https://en.wikipedia.org/wiki/Facility_location_problem): Greedy solver for the facility location problem

## Systems Analysis and Simulation - Inventory and Project Management
- [Economic Order Quantity](https://pt.wikipedia.org/wiki/Economic_order_quantity): Deterministic model featuring several variations of the original model namely 1) quantity discounts; 2) stock rupture allowed; 3) mixed stock rupture and quantity discounts.
- [Stochastic Model with Continuous Review](https://link.springer.com/chapter/10.1007/978-3-642-87146-7_4): Estimates the parameters of an economic model with stochastic demand and simulates a system's behaviour for the indicated period of time. Can be set to solve for maximum level of service or minimum cost per period.
- [Critical Path Method](https://en.wikipedia.org/wiki/Critical_path_method): Algorithm for scheduling a set of project activities. A critical path is determined by identifying the longest stretch of dependent activities and measuring the time required to complete them from start to finish. Plots a final schedule with critical and non-critical activities and respective slacks.
- [Project Evaluation and Review Technique](https://en.wikipedia.org/wiki/Program_evaluation_and_review_technique): tool used in project management to analyze and represent tasks involved in completing a given project. Estimates how much time a project will take given the duration (stochastic) of its activities. Used in conjunction with the critical path method.

## Multivariate Exploratory Analysis - Dimensionality Reduction and Clustering
- [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis): A function that provides a full detailed report on the data based on a PCA for an informed process of dimensionality reduction and data interpretation.
- Clustering Methods: From [hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering) to [K-means clustering](https://en.wikipedia.org/wiki/K-means_clustering), separates and validates the classification with optional plots for better visualization
