
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix


# 1. Read fragments from the file
with open("source.txt", "r") as file:
    numbers = [line.strip() for line in file if line.strip()]

if not numbers or any(not n.isdigit() or len(n) < 2 for n in numbers):
    raise ValueError("Invalid input data")


# 2. Create graph vertices and edges
vertices = sorted(set(
    [n[:2] for n in numbers] +
    [n[-2:] for n in numbers]
))

vertex_id = {v: i for i, v in enumerate(vertices)}

edges = [
    (vertex_id[n[:2]], vertex_id[n[-2:]])
    for n in numbers
]

m = len(edges)
n = len(vertices)


# 3. Create variables
x = np.arange(m)
start = np.arange(m, m + n)
end = np.arange(m + n, m + 2 * n)
active = np.arange(m + 2 * n, m + 3 * n)

forward = np.arange(m + 3 * n, m + 3 * n + m)
backward = np.arange(m + 3 * n + m, m + 3 * n + 2 * m)
source = np.arange(m + 3 * n + 2 * m, m + 4 * n + 2 * m)

variable_count = 3 * m + 4 * n
M = n

rows = []
lower_bounds = []
upper_bounds = []


# Helper function for adding constraints
def add_constraint(coefficients, lower=-np.inf, upper=np.inf):
    rows.append(coefficients)
    lower_bounds.append(lower)
    upper_bounds.append(upper)


# 4. Select exactly one start and one end vertex
add_constraint({int(i): 1 for i in start}, 1, 1)
add_constraint({int(i): 1 for i in end}, 1, 1)


# 5. Balance incoming and outgoing edges
for v in range(n):
    coefficients = {
        int(start[v]): -1,
        int(end[v]): 1
    }

    for e, (u, w) in enumerate(edges):
        if u == v:
            coefficients[int(x[e])] = coefficients.get(int(x[e]), 0) + 1

        if w == v:
            coefficients[int(x[e])] = coefficients.get(int(x[e]), 0) - 1

    add_constraint(coefficients, 0, 0)


# 6. Connect selected edges with active vertices
incident = [[] for _ in range(n)]

for e, (u, w) in enumerate(edges):
    incident[u].append(e)
    incident[w].append(e)

    add_constraint({
        int(x[e]): 1,
        int(active[u]): -1
    }, upper=0)

    add_constraint({
        int(x[e]): 1,
        int(active[w]): -1
    }, upper=0)

for v in range(n):
    coefficients = {int(active[v]): 1}

    for e in set(incident[v]):
        coefficients[int(x[e])] = -1

    add_constraint(coefficients, upper=0)


# 7. Add flow constraints to ensure connectivity
for e in range(m):
    add_constraint({
        int(forward[e]): 1,
        int(x[e]): -M
    }, upper=0)

    add_constraint({
        int(backward[e]): 1,
        int(x[e]): -M
    }, upper=0)

for v in range(n):
    add_constraint({
        int(source[v]): 1,
        int(start[v]): -M
    }, upper=0)

# The source sends one unit of flow to every active vertex
coefficients = {int(source[v]): 1 for v in range(n)}

for v in range(n):
    coefficients[int(active[v])] = -1

add_constraint(coefficients, 0, 0)


# Flow balance at every vertex
for v in range(n):
    coefficients = {
        int(source[v]): 1,
        int(active[v]): -1
    }

    for e, (u, w) in enumerate(edges):
        if u == v:
            coefficients[int(forward[e])] = -1
            coefficients[int(backward[e])] = 1

        if w == v:
            coefficients[int(forward[e])] = 1
            coefficients[int(backward[e])] = -1

    add_constraint(coefficients, 0, 0)


# 8. Convert constraints into a matrix
row_indices = []
column_indices = []
values = []

for row_number, row in enumerate(rows):
    for column, value in row.items():
        if value:
            row_indices.append(row_number)
            column_indices.append(column)
            values.append(value)

matrix = coo_matrix(
    (values, (row_indices, column_indices)),
    shape=(len(rows), variable_count)
).tocsr()


# 9. Define variable limits
lower = np.zeros(variable_count)
upper = np.full(variable_count, np.inf)

upper[x] = 1
upper[start] = 1
upper[end] = 1
upper[active] = 1

upper[forward] = M
upper[backward] = M
upper[source] = M

integrality = np.zeros(variable_count)
integrality[np.concatenate((x, start, end, active))] = 1


# 10. Maximize the total length of the sequence
objective = np.zeros(variable_count)

for i, number in enumerate(numbers):
    objective[i] = -(len(number) - 2)

result = milp(
    objective,
    integrality=integrality,
    bounds=Bounds(lower, upper),
    constraints=LinearConstraint(
        matrix,
        np.array(lower_bounds),
        np.array(upper_bounds)
    ),
    options={"mip_rel_gap": 0.0}
)

if not result.success:
    raise RuntimeError("Could not find an optimal solution: " + result.message)


# 11. Get selected fragments
selected = [
    i for i in range(m)
    if result.x[x[i]] > 0.5
]

start_vertex = int(np.argmax(result.x[start]))


# 12. Restore the correct order of fragments
adjacency = [[] for _ in range(n)]

for edge_id in selected:
    u, w = edges[edge_id]
    adjacency[u].append(edge_id)

vertex_stack = [start_vertex]
edge_stack = []
path = []

while vertex_stack:
    current = vertex_stack[-1]

    if adjacency[current]:
        edge_id = adjacency[current].pop()
        next_vertex = edges[edge_id][1]

        vertex_stack.append(next_vertex)
        edge_stack.append(edge_id)

    else:
        vertex_stack.pop()

        if edge_stack:
            path.append(edge_stack.pop())

path.reverse()

if len(path) != len(selected):
    raise RuntimeError("Could not restore the complete sequence")


# 13. Build the final sequence
chain = [numbers[i] for i in path]

answer = chain[0]

for number in chain[1:]:
    answer += number[2:]


# 14. Print the result
print("Total fragments:", len(numbers))
print("Used fragments:", len(chain))
print("Sequence length:", len(answer))

print("\nLongest sequence:")
print(answer)

print("\nFragment order:")
print(" -> ".join(chain))
