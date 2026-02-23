<!--Grid组件-->

<template>
    <div class="GridVis1-container">
        <svg ref="canvas" class="canvas"></svg>
        <div class="GridVis1-controller-bar">
            <el-select v-model="GeneList" multiple style="width:230px" size="mini" :collapse-tags="true">
                <el-option
                    v-for="gene in GeneList"    
                    :key="gene"
                    :label="gene"
                    :value="gene"
                ></el-option>
            </el-select>
            <el-autocomplete
                v-model="ChoosenGene"
                size="mini"
                :fetch-suggestions="geneQuerySearch"
                @select="handleGeneSelect"
                style="width:300px;">
            </el-autocomplete>
            <el-button type="primary" size="mini" @click="handleClickEnterButton">确定</el-button>
        </div>
    </div>
</template>

<script>
import * as d3 from 'd3';
import Vue from "vue";
import axios from "axios";

import {Select,Autocomplete,Button} from "element-ui";

Vue.component(Select.name, Select);
Vue.component(Autocomplete.name, Autocomplete);
Vue.component(Button.name, Button);

export default {
    name:'GridVis1',
    props:['drawData'],
    data(){
        return {
            //config
            width:800,
            height:800,
            padding:{
                top:20,
                right:20,
                bottom:20,
                left:20
            },
            patchSize:15,

            //data
            GeneList:[],
            ChoosenGene:'',


        }
    },
    methods:{
        draw(){//绘图
            const self = this;
            const svg = d3.select(this.$refs.canvas).attr('width',this.width).attr('height',this.height);
            svg.selectAll('*').remove()
            console.log(this.drawData);

            const grid_pos = self.drawData.grid_pos;
            const grid_size = self.drawData.grid_size;
            const grid_bounds = self.drawData.grid_bounds;
            const GeneVeloDirCI = self.drawData.GeneVeloDirCI;
            const meanVelo = self.drawData.meanVelo;
            const grid_colors = self.drawData.grid_colors
            const patches_data = [];
            for(let i=0;i<grid_pos.length;i++){
                patches_data.push({
                    pos:grid_pos[i],
                    dirCI:GeneVeloDirCI[i],
                    meanVelo:meanVelo[i],
                    color:grid_colors[i],
                });
            }

            //按照meanVelo的范围，确定四种大小


            console.log('patches_data:',patches_data);

            const minX = Math.min(...grid_pos.map(d=>d[0])) - 0.5 * grid_size;
            const maxX = Math.max(...grid_pos.map(d=>d[0])) + 0.5 * grid_size;
            const minY = Math.min(...grid_pos.map(d=>d[1])) - 0.5 * grid_size;
            const maxY = Math.max(...grid_pos.map(d=>d[1])) + 0.5 * grid_size;

            const xScale = d3.scaleLinear()
                .domain([grid_bounds[0][0], grid_bounds[0][1]])
                .range([this.padding.left, this.width - this.padding.right]);
            
            const yScale = d3.scaleLinear()
                .domain([grid_bounds[1][0], grid_bounds[1][1]])
                .range([this.height - this.padding.bottom, this.padding.top]);
            
            const patches_layer = svg.append('g')
            const patches = patches_layer.selectAll('g')
                .data(patches_data)
                .join('g')
                .attr("id",(d,i)=>`grid${i}`)
                .attr('transform', d=>`translate(${xScale(d['pos'][0])},${yScale(d['pos'][1])})`);

            patches.append('circle')
                .attr('r', grid_size/2 * (xScale(1)-xScale(0)))
                .attr('fill-opacity',0)
                .attr('stroke',d=>d.color)
                .attr('stroke-width',2);
            
            // arc
            const DirCIarc = d3.arc().innerRadius(0).outerRadius(grid_size/2 * (xScale(1)-xScale(0)));
            patches.selectAll('g')
                .data(d=>d.dirCI)
                .join('path')
                .attr('d', d=>DirCIarc({
                    startAngle:d['interval'][0],
                    endAngle:d['interval'][1]
                }))
                .attr('fill','black')
                .attr('opacity',d=>{
                    return d['strength'];
                })
                

            //mean velo (arrow)
            patches
                .append('line')
                .attr('x1',0)
                .attr('y1',0)
                .attr('x2', d => {
                    const vx = d.meanVelo[0];
                    const vy = d.meanVelo[1];
                    const len = Math.hypot(vx, vy) || 1e-6; // 防止除 0

                    const radiusPx =
                    0.5 * grid_size *
                    (xScale.range()[1] - xScale.range()[0]) /
                    (grid_bounds[0][1] - grid_bounds[0][0]);

                    return (vx / len) * radiusPx;
                })
                .attr('y2', d => {
                    const vx = d.meanVelo[0];
                    const vy = d.meanVelo[1];
                    const len = Math.hypot(vx, vy) || 1e-6;

                    const radiusPx =
                    0.5 * grid_size *
                    (yScale.range()[0] - yScale.range()[1]) /
                    (grid_bounds[1][1] - grid_bounds[1][0]);

                    return -(vy / len) * radiusPx;
                })
                .attr('stroke','red')
                .attr('stroke-width',1.5)
                // .attr('marker-end','url(#arrow)');
            //定义箭头
            // const defs = svg.append("defs");
            // defs.append("marker")
            //     .attr("id", "arrow")
            //     .attr("viewBox", "0 -5 10 10")
            //     .attr("refX", 10)          // 箭头尖端对齐线的终点
            //     .attr("refY", 0)
            //     .attr("markerWidth", 6)
            //     .attr("markerHeight", 6)
            //     .attr("orient", "auto")
            //     .append("path")
            //     .attr("d", "M0,-5L10,0L0,5")
            //     .attr("fill", "red");

        },
        geneUpdatePlot(){//按照基因，绘图更新

        },
        geneQuerySearch(queryString, cb){

            let matched_genes = this.drawData.Genes.filter(str => str.toLowerCase().startsWith(queryString.toLowerCase()))
            let sorted_matched_genes = matched_genes.sort((a,b)=>a.length - b.length)
            cb(sorted_matched_genes.map(item=>{
                return { value: item };
            }))
        },
        handleGeneSelect(item){
            if(this.GeneList.includes(item.value)){//基因已经有了该基因
                this.$message({
                    'message':'This Gene has been added',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            this.GeneList.push(item.value);
        },
        handleClickEnterButton(){
            axios({
                    method:"post",
                    url:"/api/update_plot_GridVis1",
                    data:{
                        'project_name':'Pancreatic_byGenes_20neighbors_2top_velo',
                        'gene_list':this.GeneList,
                    }
                }).then(res=>{
                    let data = res.data;
                    // 更新GeneVeloDirCI和meanVelo
                    this.drawData.GeneVeloDirCI = data.GeneVeloDirCI;
                    this.drawData.meanVelo = data.meanVelo;
                    
                    this.draw();
                    console.log('data:',data)
                }).catch(err=>{
                    console.log('err:',err)
                })
                this.geneUpdatePlot();

        }
    },
    watch:{
        drawData:{
            handler(){
                this.draw();
            },
            deep:true,
        },
    },
    mounted(){
        this.draw();
    }
    
}
</script>

<style>
.GridVis1-container{
    display: flex;
    flex-direction: column;

}

.controller-bar{
    width:100%;
    height:40px;
    border-bottom:1px solid #ccc;
    box-sizing:border-box;
}

</style>