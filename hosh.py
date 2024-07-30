'''
Artificial intelligence project ~ solve the Knapsack by Genetic Algorithm
Fatemeh Rahmani ~ 990201110022'''

# Import libraries
from typing import List
from random import random
from colorama import Fore, Style
import copy
import csv

# Name the objects inside the backpack.
products_list = open('products_list.csv', 'r') 
reader = csv.reader(products_list) 
list = []  
for record in reader: 
    list.append(record)  

# We consider the values ​​of each line as a list
name = [float(i) for i in list[0][1:]]
space = [float(i) for i in list[1][1:]] 
price = [float(i) for i in list[2][1:]]

# When we want to create a class later and so
#  that the interpreter doesn't give an error, we put this value in it
class backpack:
    pass 


class Product:  
    # It is similar to a gene in the genetic algorithm.
    def __init__(self, name:str, space:float, price:int) -> None:
        self.name = name
        self.space = space
        self.price = price
        
    def __repr__(self) -> str:
        return f"{self.name} * {self.space} * {self.price}"


# Specifying a class for genetic algorithm chromosomes and defining them as well as updating them by periods
class backpack:# It is similar to the chromosome in the genetic algorithm.
    def __init__(self, products:List[Product], space_limit:float, generation:int=0) -> None:
        self.products = products                # The maximum weight it can bear
        self.space_limit = space_limit
        self.generation = generation
        self.chromosome = self.made_chromosome() 
        self._selected_products = None
        self.reperantion = None
        self.evaluated_score = 0 # Price
        self.used_space = 0
        
        self.update_backpack_info()
        
    def made_chromosome(self) -> List[int]:
        chromosome = []
        for _ in range(len(self.products)):
            if random() > 0.5:
                chromosome.append(0) # Add 0
            else:
                chromosome.append(1) # Add 1    
        return chromosome
    
    def create_one_child(self, list_new_chromose:List[int]) -> backpack:
       # Creates the next generation child with the new chromosome.
        child_our_parents = backpack(self.products, self.space_limit, self.generation + 1)
        child_our_parents.chromosome = list_new_chromose
        child_our_parents.update_backpack_info()
        
        return child_our_parents
        
    def __gt__(self, _other_backpack:backpack) -> bool:
        return self.evaluated_score > _other_backpack.evaluated_score
    
    def __eq__(self, _other_backpack:backpack) -> bool:
        return self.selected_products == _other_backpack.selected_products

    def __repr__(self) -> str:
        # Total reperantion 
        return self.reperantion    

    @property  # The function used to get the attribute value
    def selected_products(self) -> List[Product]:
        return self._selected_products

    def selected_products_update(self) -> None:
        selected_products_list = []
        for gene, product  in zip(self.chromosome, self.products): # Matches chromosomes and products pairwiseZIP
            if gene == 1:
                selected_products_list.append(product)  # Adds the product to the list of selected products
                
        self._selected_products = selected_products_list
        
    def Structure_of_output_answer(self) -> None:
        string_reperantion = ""
        seperator = "======================================================\n" # The dividing line
        string_reperantion += seperator 

        for idx, product in enumerate(self.selected_products): # enumerate iterates through the list
            product_desc = f"Product {idx+1}: {repr(product)} \n" # Product Name
            string_reperantion += product_desc
        
        string_reperantion += f"Generation ~> {self.generation} "
        string_reperantion += f"And the Chromosome ~> {self.chromosome} \n"
        string_reperantion += f"Score ~> {self.evaluated_score} \n"

        string_reperantion += seperator
        self.reperantion = string_reperantion
    
    def update_total_fitness_info(self) -> None:
        #We consider weight and price as fitness.
        total_space = 0
        total_score = 0 # price
        for product in self.selected_products:
            total_space += product.space
            total_score += product.price 
            
        if self.space_limit < total_space:
            total_score = 0
        self.used_space = total_space
        self.evaluated_score = total_score
        
    def update_backpack_info(self) -> None:
      # Display information in the output
        self.selected_products_update()
        self.update_total_fitness_info()
        self.Structure_of_output_answer()
    
    def crossover(self, another_backpack:backpack) -> List[backpack]:
        number_of_cut_off_point = n
        # math.floor(random() * len(self.chromosome)) # To maintain diversity

        self_genes_x, self_genes_y = self.chromosome[:number_of_cut_off_point], self.chromosome[number_of_cut_off_point:]
        other_genes_x, other_genes_y = another_backpack.chromosome[:number_of_cut_off_point],  another_backpack.chromosome[number_of_cut_off_point:]  
        
        chromosome_child_x = [*self_genes_x, *other_genes_y]
        chromosome_child_y = [*other_genes_x, *self_genes_y]  
        
        next_generation_child_x = self.create_one_child(chromosome_child_x)
        next_generation_child_y = self.create_one_child(chromosome_child_y)
        
        children_of_next_generation = [next_generation_child_x, next_generation_child_y]
        return children_of_next_generation
        
    #Mutation of a gene in the chromosome
    def mutate_one_gene(self, idx:int) -> None:
        bit = self.chromosome[idx]
        newly_built_bit = 0 if bit == 1 else 1
        self.chromosome[idx] = newly_built_bit
        
    # Chromosome mutation
    def mutation(self, mutation_rate:float) -> backpack:
        #print(f"Before: {self.chromosome}")
        mutated = False
        for idx in range(len(self.chromosome)):
            mutation_prob = random()
            if mutation_prob > mutation_rate:
                continue
                
            self.mutate_one_gene(idx)
            mutated = True

        if mutated:
            self.update_backpack_info()
            #print(f"After: {self.chromosome}")
        return self

# made Without using the library
class Genetic_Algorithm():
    def __init__(self, population_size:int, mutation_prob, num_of_generation, products:List[Product], space_limit:int) -> None:
        self.population_size = population_size
        self.population = [].copy()
        self.generation = 0
        self.number_of_generation = num_of_generation
        self.find_best_solution = None
        self.find_global_best = None
        self.mutation_probability = mutation_prob
        self.products = products
        self.space_limit = space_limit
    
    def show_header(self) -> None:
        print(Fore.LIGHTCYAN_EX ,f"checking ...\n", Fore.RESET, f"  Generation: {self.find_best_solution.generation}",f"And the Best Solution is: \n{self.find_best_solution}")
        
    # Initial population
    def primary_population(self) -> None:
        for _ in range(self.population_size):
            one_solution = backpack(self.products, self.space_limit)
            self.population.append(one_solution)
        self.population.sort(reverse = True)
        self.find_best_solution = self.population[0]
        self.find_global_best = self.find_best_solution
        
    def t_fitness_population(self) -> int:
       #Our people come and calculate the value of fitness and add it to the roulette wheel.
        t_fitness = 0
        for backpack in self.population:
            t_fitness += backpack.evaluated_score    
        return t_fitness
    
    # Auxiliary function for making the roulette wheel
    def select_partent(self) -> backpack:
        t_fitness = self.t_fitness_population()
        random_value = random() * t_fitness
        
        check_summation = 0
        what_selected_idx = 0
        #print(f"Random Value: {random_value}")

        while (self.population_size > what_selected_idx and random_value > check_summation):
            #print(f"Selected Chromosome Index: {selected_idx}, Check Sum: {check_sum}")
            check_summation += self.population[what_selected_idx].evaluated_score # ارزیابی
            what_selected_idx += 1
            #print(f"Selected Chromosome Index: {selected_idx}, Check Sum: {check_sum}")
        return self.population[what_selected_idx - 1]

    def crossover_parents(self) -> List[backpack]:
        parent_a = self.select_partent()
        parent_b = self.select_partent()
        while (parent_a == parent_b):
            parent_b = self.select_partent() # Selects another parent.
        make_children = parent_a.crossover(parent_b)
        return make_children
    
    def mutation_children(self,mutation_prob:float,children:List[backpack]) -> List[backpack]:
        copy_children = copy.deepcopy(children)
        for children in copy_children:
            children.mutation(mutation_prob)    
        return copy_children
        
    
    def solve_algorithm(self) -> backpack:
        global generation
        self.primary_population() 
        self.show_header()
        
        for generation in range(self.number_of_generation):
            new_population = []
            for _ in range(self.population_size // 2): # To create a parent.
                # crossover
                child_x, child_y = self.crossover_parents()
                # mutation
                mutated_child_x, mutated_child_y = self.mutation_children(self.mutation_probability, [child_x, child_y])

                new_population.extend([mutated_child_x, mutated_child_y])
            self.population = new_population
            self.population.sort(reverse = True)
            self.find_best_solution = self.population[0]
            self.find_global_best = max(self.find_global_best, self.find_best_solution)
            self.show_header()   
        return self.find_global_best


def main():
    global n
    print(Fore.GREEN, Style.BRIGHT,"Hello :)\nThank you for your choice.")
    print("This program solves the backpack problem by genetic algorithm and displays partially optimal answers.\n", Style.NORMAL, Fore.RESET)
    
    # Get input
    population_size = int(input("Please enter the population size (1 to 14) you want:"))
    mutation_prob = float(input("Please enter the mutation probability you want:"))
    num_of_generation = int(input("Please enter the num of generation you want:"))
    products_l = products_list
    # Maximum weight of the bag
    space_limit = int(input("Please enter the maximum weight you want your backpack to accept:"))
    elitism = int(input("if you want apply elitism press 1 If not, press 0:"))
    n = int(input("Please enter the number of 'n' for n point crossover:")) 
    print("\n")

    genetic_algorithm = Genetic_Algorithm(population_size, mutation_prob, num_of_generation, products_l, space_limit)
    genetic_algorithm.solve_algorithm()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('ERROR')
      