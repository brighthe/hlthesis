此外，我发现 Marp 2. 可以逐步显示内容
你可以使用 CSS 动画来逐步显示幻灯片内容：

```markdown
---
marp: true
style: custom-style.css
---

# 逐步显示内容示例

<div class="step-by-step">
    <p class="step">第一步</p>
    <p class="step">第二步</p>
    <p class="step">第三步</p>
</div>

<style>
@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

.step-by-step .step {
    opacity: 0;
    animation: fadeIn 1s forwards;
    animation-delay: calc(1s * var(--step));
}

.step-by-step .step:nth-child(1) { --step: 0; }
.step-by-step .step:nth-child(2) { --step: 1; }
.step-by-step .step:nth-child(3) { --step: 2; }
</style>
```