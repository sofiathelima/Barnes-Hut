package main

import (
	"math"
)

func (q *Quadrant) NewNode() *Node {

	n := &Node{
		children: make([]*Node, 4),
		sector:   *q,
	}

	subquads := q.Partition()

	for i := range n.children {
		n.children[i] = &Node{
			sector: subquads[i],
		}
	}
	return n
}

func (n *Node) SetStar(leaf *Star) *Node {
	n.star = leaf
	return n
}

func NewQuadrant(_x, _y, w float64) Quadrant {
	return Quadrant{
		x:     _x,
		y:     _y,
		width: w,
	}
}

// Partition takes a Quadrant, and returns a list of 4 subquadrants
// NE, NW, SE, SW
func (quad *Quadrant) Partition() []Quadrant {
	_x := quad.x
	_y := quad.y
	w := quad.width

	NE := NewQuadrant(_x+w/2, _y+w/2, w/2)
	NW := NewQuadrant(_x, _y+w/2, w/2)
	SE := NewQuadrant(_x+w/2, _y, w/2)
	SW := NewQuadrant(_x, _y, w/2)

	return []Quadrant{NE, NW, SE, SW}
}

// SumMasses ranges a list of stars and accumulates mass
func SumMasses(s []*Star) float64 {
	var sum float64
	for _, star := range s {
		if star != nil {
			sum += star.mass
		}
	}
	return sum
}

// SetWeightedPosition takes a list of starts and calculates a
// weighted center of mass position.
func SetWeightedPosition(s []*Star) OrderedPair {

	var _x float64
	var _y float64
	for _, star := range s {
		if star != nil {
			_x += star.position.x * star.mass
			_y += star.position.y * star.mass
		}
	}
	masses := SumMasses(s)

	return OrderedPair{x: _x / masses, y: _y / masses}
}

// SetDummyStar acts on a *Node. It creates a list of stars
// from the stars of the Node's 4 children nodes. Then uses
// that list to calculate
func (n *Node) SetDummyStar() *Node {

	childStars := make([]*Star, 0)
	for _, child := range n.children {
		if child.star != nil {
			childStars = append(childStars, child.star)
		}
	}

	dummy := Star{}
	dummy.mass = SumMasses(childStars)
	dummy.position = SetWeightedPosition(childStars)
	n.SetStar(&dummy)

	return n
}

func (n *Node) MapStarsToSubQuads(stars []*Star) map[Quadrant][]*Star {

	// values of current Quadrant
	_x := n.sector.x
	_y := n.sector.y
	w := n.sector.width

	// Make a map of subquadrants(stored in node's children) to empty list of stars

	NE_stars := make([]*Star, 0)
	NW_stars := make([]*Star, 0)
	SE_stars := make([]*Star, 0)
	SW_stars := make([]*Star, 0)

	// Range all stars and append list
	for _, star := range stars {
		if star.position.y > _y+w/2 { // if in top half
			if star.position.x > _x+w/2 {
				NE_stars = append(NE_stars, star) // children in top right (NE)
			} else {
				NW_stars = append(NW_stars, star) // NW
			}
		} else {
			if star.position.x > _x+w/2 {
				SE_stars = append(SE_stars, star) // children in bottom right (SE)
			} else {
				SW_stars = append(SW_stars, star) // SW
			}
		}

	}

	dividedStars := map[Quadrant][]*Star{
		n.children[0].sector: NE_stars,
		n.children[1].sector: NW_stars,
		n.children[2].sector: SE_stars,
		n.children[3].sector: SW_stars,
	}

	return dividedStars
}

func (n *Node) RecursiveBuildQuadTree2(stars []*Star, currQuad Quadrant) *Node {

	if len(stars) == 0 {
		return n.SetStar(nil)
	}
	if len(stars) == 1 {
		return n.SetStar(stars[0])
	}
	newNode := currQuad.NewNode()
	mappedStars := newNode.MapStarsToSubQuads(stars)

	// range through children of node, use the list of stars mapped to the sector
	// of each child, and return the node of this quadrant
	call := 0
	for i, child := range newNode.children {
		call++
		newNode.children[i] = child.RecursiveBuildQuadTree2(mappedStars[child.sector], child.sector)
	}

	return newNode.SetDummyStar()

}

// NewPosition computes the new poosition given the updated acc and velocity.
//
// Assumputions: constant acceleration over a time step.
// => DeltaX = v_avg * t
//    DeltaX = (v_start + v_final)*t/ 2
// because v_final = v_start + acc*t:
//	  DeltaX = (v_start + v_start + acc*t)t/2
// Simplify:
//	DeltaX = v_start*t + 0.5acc*t*t
// =>
//  NewX = v_start*t + 0.5acc*t*t + OldX
//
func (b *Star) NewPosition(t float64) OrderedPair {
	return OrderedPair{
		x: b.position.x + b.velocity.x*t + 0.5*b.acceleration.x*t*t,
		y: b.position.y + b.velocity.y*t + 0.5*b.acceleration.y*t*t,
	}
}

// Compute the Euclidian Distance between two bodies
func Dist(b1, b2 *Star) float64 {
	dx := b1.position.x - b2.position.x
	dy := b1.position.y - b2.position.y
	return math.Sqrt(dx*dx + dy*dy)
}

// ComputeGravityForce computes the gravity force between body 1 and body 2.
func ComputeGravityForce(b1, b2 *Star) OrderedPair {
	d := Dist(b1, b2)
	deltaX := b2.position.x - b1.position.x
	deltaY := b2.position.y - b1.position.y
	F := G * b1.mass * b2.mass / (d * d)

	return OrderedPair{
		x: F * deltaX / d,
		y: F * deltaY / d,
	}
}

func (v *OrderedPair) Add(v2 OrderedPair) {
	v.x += v2.x
	v.y += v2.y
}

// ComputeNetForce sums the forces of all bodies in the universe
// acting on b.
func ComputeNetForce(univ *Universe, b *Star, theta float64, tree QuadTree) OrderedPair {
	var netForce OrderedPair
	for _, body := range univ.stars {
		if body != b {
			f := tree.root.TreeForce(body, theta)
			netForce.Add(f)
		}
	}
	return netForce
}

func (n *Node) TreeForce(body *Star, theta float64) OrderedPair {
	var netForce OrderedPair
	s := n.sector.width

	// leaf --> compute foce
	// internal node --> check s/d8u

	for _, child := range n.children {
		if child.star == body {
			break
		}
		if child != nil {
			if s/Dist(body, n.star) < theta {
				return ComputeGravityForce(body, n.star)
			} else {
				netForce.Add(child.TreeForce(child.star, theta))
				return netForce
			}
		}
	}

	return netForce
}

// UpdateAccel computes the new accerlation vector for b
func (b *Star) NewAccel(univ *Universe, theta float64, tree QuadTree) OrderedPair {
	F := ComputeNetForce(univ, b, theta, tree)
	return OrderedPair{
		x: F.x / b.mass,
		y: F.y / b.mass,
	}
}

// NewVelocity makes the velocity of this object consistent with the acceleration.
func (b *Star) NewVelocity(t float64) OrderedPair {
	return OrderedPair{
		x: b.velocity.x + b.acceleration.x*t,
		y: b.velocity.y + b.acceleration.y*t,
	}
}

// need to change this to use current acc/vel/pos
func (b *Star) Update(univ *Universe, t, theta float64, tree QuadTree) {

	acc := b.NewAccel(univ, theta, tree)
	vel := b.NewVelocity(t)
	pos := b.NewPosition(t)
	b.acceleration, b.velocity, b.position = acc, vel, pos
}

// UpdateUniverse returns a new universe after time t.
func UpdateUniverse(univ *Universe, t, theta float64, tree QuadTree) *Universe {

	newUniverse := CopyUniverse(univ)

	for b := range univ.stars {

		// update pos, vel and accel
		newUniverse.stars[b].Update(univ, t, theta, tree)
	}

	return newUniverse
}

func (g *Galaxy) Push(x, y float64) {
	g.g[len(g.g)-1].velocity.x = x
	g.g[len(g.g)-1].velocity.x = y
}

//BarnesHut is our highest level function.
//Input: initial Universe object, a number of generations, and a time interval.
//Output: collection of Universe objects corresponding to updating the system
//over indicated number of generations every given time interval.
func BarnesHut(initialUniverse *Universe, numGens int, time, theta float64) []*Universe {
	timePoints := make([]*Universe, numGens+1)
	timePoints[0] = initialUniverse

	for gen := 0; gen < numGens; gen++ {
		var tree QuadTree
		quad := NewQuadrant(0, 0, timePoints[gen].width)
		tree.root = tree.root.RecursiveBuildQuadTree2(timePoints[gen].stars, quad)

		timePoints[gen+1] = UpdateUniverse(timePoints[gen], time, theta, tree)

	}

	return timePoints
}
