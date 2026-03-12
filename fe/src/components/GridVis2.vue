<!--简单的embedding + 二维速度箭头的组件-->

<template>
    <div class="GridVis2-container">
        <svg ref="canvas" class="canvas"></svg>
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
    name:'GridVis2',
    props:['drawData'],
    data(){
        return {
            //config
            width:800,
            height:800,
            scatter_size:3,
            padding:{
                top:20,
                right:20,
                bottom:20,
                left:20
            },
            maxArrowLength:35,
        }
    },
    methods:{
        draw(){//绘图
            const self = this;
            const svg = d3.select(this.$refs.canvas).attr('width',this.width).attr('height',this.height);
            svg.selectAll('*').remove()
            console.log(this.drawData);

            //整理数据
            let scatter_data = []
            for(let i = 0; i < this.drawData['embedding'].length; i++){
                scatter_data.push({
                    'pos': this.drawData['embedding'][i],
                    'velocity': this.drawData['velocity_embedding'][i],
                    'color':this.drawData['cluster_color'][i],
                })
            }
            //按照meanVelo的范围，确定四种大小
            console.log('scatter_data:',scatter_data);

            const minX = Math.min(...scatter_data.map(d=>d['pos'][0]));
            const maxX = Math.max(...scatter_data.map(d=>d['pos'][0]));
            const minY = Math.min(...scatter_data.map(d=>d['pos'][1]));
            const maxY = Math.max(...scatter_data.map(d=>d['pos'][1]));
            const maxVelocity = d3.max(scatter_data, d =>
                Math.hypot(d.velocity[0], d.velocity[1])
            ) || 1e-6;

            const xScale = d3.scaleLinear()
                .domain([minX, maxX])
                .range([this.padding.left, this.width - this.padding.right]);
            
            const yScale = d3.scaleLinear()
                .domain([minY, maxY])
                .range([this.height - this.padding.bottom, this.padding.top]);
            
            const scatter_layer = svg.append('g')

            const scatter = scatter_layer.selectAll('g')
                .data(scatter_data)
                .join('g')
                .attr('transform', d=>`translate(${xScale(d['pos'][0])},${yScale(d['pos'][1])})`);
            
            //dot
            const dot = scatter.append('circle')
                .attr('r', this.scatter_size)
                .attr('fill',d=>d.color)

            //arrow
            const arrowScale = d3.scaleSqrt()
                .domain([0, maxVelocity])
                .range([0, 1]) // 技巧：将 range 的起点设为 0.2 而不是 0，保证最小速率也有个可见基数
                .clamp(true);
            const arrow = scatter.append('line')
                .attr('x1', 0)
                .attr('y1', 0)
                .attr('x2', d => {
                    const vx = d.velocity[0];
                    const vy = d.velocity[1];

                    const len = Math.hypot(vx, vy) || 1e-6;

                    // 归一化到 [0, 3]

                    const scale = arrowScale(len)
                    return self.maxArrowLength * scale * (vx / len);
                })
                .attr('y2', d => {
                    const vx = d.velocity[0];
                    const vy = d.velocity[1];

                    const len = Math.hypot(vx, vy) || 1e-6;
                    
                    const scale = arrowScale(len)
                    
                    
                    return -self.maxArrowLength * scale * (vy / len);
                })
                .attr('stroke', 'black')
                .attr('stroke-width', 1.5);

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
.GridVis2-container{
    display: flex;
    flex-direction: column;

}

</style>