package main

import . "github.com/scipipe/scipipe"

func main2() {
	wf := NewWorkflow("hello_workflow", 4)

	hello := wf.NewProc("hello", "echo 'Hello, world!' > {o:out}")
	hello.SetOut("out", "hello.txt")

	wf.Run()
}
