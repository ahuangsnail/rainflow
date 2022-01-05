namespace RainFlow

type binsInpType =
    | Bin of nbin:int
    | BinsWithBoundary of xmin:double * xmax:double * nbin:int

type HistogramSeriesType = 
    | Xs of xs:array<double>
    | XsByWeight of xs:array<double>*weights:array<double>

[<AutoOpen>]
module RainFlow =
    type IndexPare = int * double
    type CycleCountedType = double * double * double * int * int

    let enumerate (seris:list<double>) (ind: int) =
        match seris with
        | [] -> []
        | _ ->
            List.mapi (fun i x -> (i+ind, x)) seris


    let reversals (serisinp) : list<IndexPare> =
        let seris = enumerate serisinp 0
        
        /// 返回seris序列的相反数和序列值
        match ((Seq.isEmpty seris)
               || (Seq.tail seris |> Seq.isEmpty)
               || (Seq.tail seris |> Seq.tail |> Seq.isEmpty)) with
        | true -> seris |> Seq.toList
        | _ ->
            let firstp: IndexPare = Seq.head seris
            let (ixFirst, xFirst) = firstp
            let (ixSecond, xSecond) = Seq.tail seris |> Seq.head
            let mutable xPreDiff = xFirst - xSecond
            let mutable ixnext = ixFirst
            let mutable xnext = xFirst
            
            let pareseq = List.windowed 2 seris
            let mutable pickSeris = [(ixFirst,xFirst)]
            for item in pareseq do
                let (ixc, xc) = item.[0]
                ixnext <- fst item.[1]
                xnext <- snd item.[1]
                let xDiff = xc - xnext

                if xPreDiff * xDiff < 0. then
                    pickSeris <- pickSeris@[(ixc, xc)]

                xPreDiff <- xDiff
            pickSeris <- pickSeris@[(ixnext, xnext)]
            pickSeris
    
    


    let countCycle (seris) :list<CycleCountedType> =
        let formatOutput point1 point2 count =
            let (i1, x1) = point1
            let (i2, x2) = point2
            let rng = abs (x1 - x2)
            let mean = 0.5 * (x1 + x2)
            rng, mean, count, i1, i2

        let mutable rfhist = []
        let mutable pointstack: array<int * double> = [| |]
        let mutable ifcontinue = true
        for pnt in reversals seris do
            pointstack <- Array.append pointstack  [|pnt |]
            ifcontinue <- true
            while (ifcontinue && (Array.length pointstack) >= 3) do
                let (ix1, x1) = pointstack.[^2]
                let (ix2, x2) = pointstack.[^1]
                let (ix3, x3) = pointstack.[^0]
                let X = abs (x3 - x2)
                let Y = abs (x2 - x1)

                if X < Y then
                    ifcontinue <- false
                elif (Array.length pointstack) = 3 then
                    rfhist <- rfhist @ [(formatOutput pointstack.[0] pointstack.[1] 0.5)]
                    pointstack <- Array.tail pointstack
                else
                    rfhist <- rfhist @ [ formatOutput pointstack.[^2] pointstack.[^1] 1.0]
                    pointstack <- Array.append pointstack.[..^3]  pointstack.[^0..]

        while Array.length pointstack > 1 do
            rfhist <- rfhist @ [formatOutput pointstack.[0] pointstack.[1] 0.5]
            pointstack <- Array.tail pointstack
        
        rfhist     

    let histogram (data: HistogramSeriesType) (bulk: binsInpType) = 
        let histBy xmin xmax nbin (ss: HistogramSeriesType) = 
            let fnorm = double(nbin)/(xmax-xmin)
            let cnts : array<double> = Array.zeroCreate nbin
            match ss with
            | Xs(xs) -> 
                xs
                |> Array.iter (fun x ->
                    if x>=xmin && x<xmax then
                        let ind = int((x-xmin)*fnorm)
                        cnts.[ind] <- cnts.[ind] + 1.0)
                cnts
            | XsByWeight(xs, weights) ->
                (xs, weights)
                ||> Array.zip
                |> Array.iter (fun (x,w)->
                        if x>=xmin && x<xmax then
                            let ind = int((x-xmin)*fnorm)
                            cnts.[ind] <- cnts.[ind]+w)
                cnts
                    
        match bulk with
        | Bin(nbin : int) ->
            let xs = match data with
                        | Xs(xs) -> xs
                        | XsByWeight(xs,weights) -> xs
            let xmin:double = Array.min xs
            let xmax:double = (Array.max xs) + 1e-16
            histBy xmin xmax nbin data

        | BinsWithBoundary(xmin:double, xmax:double, nbin: int) ->
            histBy xmin xmax nbin data
