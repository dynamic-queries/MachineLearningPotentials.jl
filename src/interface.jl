abstract type AbstractDescriptor end 
abstract type AbstractLearner end 

struct POD <: AbstractDescriptor end 
struct LuOpt <: AbstractLearner end 


struct MLPotentialProblem 
    trajectories::Array
    velocities::Array
end 

mutable struct MLPotentialSolution
    potential
    descriptor::AbstractDescriptor
    learner::AbstractLearner
end

function __solve(prob::MLPotentialProblem, sol::MLPotentialSolution{Any,POD,LuOpt})
    @unpack trajectories, velocities = prob
    @unpack potential = sol
end 

function solve(prob::MLPotentialProblem, desciptor::AbstractDescriptor, learner::AbstractLearner)
    solution =  MLPotentialSolution(nothing,desciptor,learner)
    __solve(prob,solution)
end